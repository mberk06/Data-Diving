#' Author: Michael Berk
#' Date: Summer 2020
#' Descrtiption: unsupervised clssification of reefs
#' 
#'     In this file we explore unsupervised classifications of Caribbean reefs using K-means.
#'     These classification models are fit using substrate values collected by Reef Check. See 
#'     post for documentation of data and collection methods. 
#'     
#'     Sections
#'     1. Helper: helper functions called later
#'     2. Exploration: unorganized exploration that informed later analysis
#'     3. Clustering Exploration: exploration centered around k-means clustering
#'     4. Final: code that was used to develop the graphs in the post
#'     
#'     Notes:
#'     - The motivation for developing clustering algorithms is the lack of a consistant TS. For dive
#'       sites there are rarely repeat dives, so we are currently aggregating using a time frame for 
#'       all reefs in the Caribbean. However, this averages out a ton of signal. If we can instead cluster
#'       all reefs by "type," we can create a more consistant time series.
#'     - The results were interesting but did not address the above question. More subject matter calls
#'       will probably need to be made to create intuitive groupings of reefs that allows for TS analysis.

pacman::p_load(ggplot2) 
pacman::p_load(reshape2) # melt
pacman::p_load(stats) # KNN
pacman::p_load(caret) # train/test
pacman::p_load(stringr) # str_replace
pacman::p_load(comprehenr) # list comps
pacman::p_load(flexclust) # weighted kmeans
pacman::p_load(factoextra) # kmeans plot
pacman::p_load(ggpubr) # multiple plots one window

# set wd based on user
if (dir.exists("/Users/michaelberk")) {
  setwd('~/Documents/Ocean Analysis/')
} else {
  setwd('~ISA PATH') # TODO CHANGE TO CORRECT PATH
}

# setup data
source(paste0(getwd(), '/Scripts/setupData.R'))

############
############
# Helpers
############
############
predict.kmeans <- function(object, newdata, isWeighted){
  #' predict cluster from trained kmeans model
  #' 
  #' @param obeject (kmeans object)
  #' @param newdata (data.frame) with new data
  #' @param isWeighted (bool) check if is weighted
  #' @return vector of clusters
  
  if (isWeighted) {
    centers <- object@centers 
  } else {
    centers <- object$centers
  }
  n_centers <- nrow(centers)
  dist_mat <- as.matrix(dist(rbind(centers, newdata)))
  dist_mat <- dist_mat[-seq(n_centers), seq(n_centers)]
  max.col(-dist_mat)
}

getSummaryStats <- function(df, cols) {
  #' calcualte summary stats for each column and return an aggregated df
  #' 
  #' @param df data.frame of all data. This should already be subset to region
  #' @param cols str vector with column names
  #' @return data.frame where columns are summary stats values and rows are cols
  
  returnDF <- data.frame(NAME=NA, Q1=NA, Q2=NA, MEAN=NA, Q3=NA)
  
  # iterte trough cols
  for (c in cols) {
    # calculate summary stats
    stats <- summary(df[,c])
    
    # append to returnDF
    returnDF <- rbind(returnDF, c(NA, stats[[2]], stats[[3]],stats[[4]], stats[[5]]))
    returnDF$NAME[dim(returnDF)[1]] <- c
  }
  
  # reset index
  rownames(returnDF) <- 1:length(returnDF[,1])
  
  return(returnDF[2:dim(returnDF)[1],])
}

substrateBarplot <- function(df, y, title, isKMeans = F, showXLabs = T, showYLabs = T) {
  #' Show stylyzed barplot with substrate names
  #' 
  #' @param df (data.frame) with data - often from allSummaries
  #' @param y (character) which column to plot
  #' @param title (character) title for plot
  #' @param isKMeans (bool) if plotting kmeans output
  #' @param showXLabs (bool) if should show xlabels
  #' @param showYLabs (bool) if should show xlabels
   
  # transform data if needed
  if (isKMeans) {
    df <- as.data.frame(t(df))
    df$NAME <- as.vector(rownames(df))
    names(df) <- c(y, 'NAME')
  }
  
  # scale y by 100
  df[,y] <- df[,y]*100
  
  # create x-axis title vector (in df$NAME)
  if (showXLabs) {
    namesList <- subsetCols('substrateNames')
    names(namesList) <- subsetCols('substrate')
    xAxisLabs <- to_vec(for(n in df$NAME) namesList[[n]])
    df$NAME <- str_replace(xAxisLabs, '\\.',' ')
  } 
  
  # create plot with fills
  p <- ggplot(data = df, aes(x=NAME, y=get(y), fill = get(y) < 0)) +
    geom_bar(stat = 'identity') + theme_minimal() +
    scale_fill_manual(guide = FALSE,
                      name = 'change < 0', 
                      values = setNames(c(red, blue), c(T, F)))
  
  # rename cols and add title
  p <- p + ggtitle(title) + 
    theme(axis.text.x = element_text(angle = 0, size = 9, hjust = 0.5), # DO 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size = 14)) 
  
  # add ylab
  if ((isKMeans) && (showYLabs)) {
    p <- p + ylab('Percent Cover Relative To Average') + xlab('')
  } else if (showYLabs) {
    p <- p + ylab('Percent Cover') + xlab('')
  } else {
    p <- p + ylab('') + xlab('')
  }
  
  # remove x lab if needed
  if (!showXLabs) {
    p <- p + theme(axis.text.x = element_blank())
  }
  
  p  
}

runKMeans <- function(df, cols, k, doTrainTest, doWeight = T, splitDate = NA) {
  #' Perform k means clustering and return percent change relative to average for each col
  #' 
  #' @param df (data.frame) data to fit
  #' @param cols (vect character) cols to cluster on
  #' @param k (numeric) integer number of clusters
  #' @param doTrainTest (bool) should perform train test split
  #' @param doWeight (bool) weights for the data
  #' @param splitDate (date) date to split on (only if doTrainTest)
  #' @return (data.frame) percent change dataframe
  
  if (doTrainTest) {
    # create train test split
    train_ind <- sample(seq_len(nrow(df)), size = floor(0.75 * nrow(df)))
    train <- na.omit(df[train_ind, cols])
    test <- df[-train_ind,]
    
    # get new and old subset for testing
    new <- na.omit(test[as.Date(test$DATE) > splitDate, cols])
    old <- na.omit(test[as.Date(test$DATE) <= splitDate, cols])
    
  } else {
    # create train (without test) for fitting
    train <- na.omit(df[,cols])
  }
  
  # calculate summary stats
  allSummaries <- getSummaryStats(train, cols)
  rownames(allSummaries) <- 1:length(allSummaries[,1])
  
  # fit kmeans with or without weights
  if (!doWeight) {
    fit <- kmeans(train, centers=k, nstart = 25)
    percentChanges <- fit$centers
  } else {
    fit <- cclust(train, k=k, weights=allSummaries$Q2, method="hardcl")
    percentChanges <- fit@centers
  }
  
  if (doTrainTest) {
    # predict based on data 
    print(table(predict.kmeans(fit, old, isWeighted = doWeight))/dim(old)[1])
    print(table(predict.kmeans(fit, new, isWeighted = doWeight))/dim(new)[1])
  } 
  
  # get percent changes from mean
  for (i in 1:dim(percentChanges)[1]) {
    for (j in 1:length(cols)) {
      meanIndex <- match(cols[j],allSummaries$NAME)
      percentChanges[i,j] <- 100*(percentChanges[i,j]-allSummaries$MEAN[meanIndex])/allSummaries$MEAN[meanIndex]
    }
  }

  return(as.data.frame(percentChanges))
}



###############
###############
# Exploration
###############
###############

# get mean, median, 25th, and 75th percentile for reefs
allSummaries <- getSummaryStats(df, substrateNames)
allSummaries

# specific value
substrateBarplot(allSummaries, 'Q2', title = 'Percent Cover for Caribbean Reefs (Median)')

# melt and plot 
summaryDF <- melt(allSummaries, id.vars='NAME')

ggplot(summaryDF, aes(x=NAME, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + 
 scale_fill_brewer()

# boxplot
boxplot(df[,subsetCols('substrate')])

# check by differnt years
allSummaries <- getSummaryStats(subset(df, (df$YEAR > 2000) & df$YEAR < 2005), substrateNames)
substrateBarplot(allSummaries, 'Q2', title = 'Percent Cover for Caribbean Reefs (Median)')

allSummaries <- getSummaryStats(subset(df, (df$YEAR > 2005) & df$YEAR < 2010), substrateNames)
substrateBarplot(allSummaries, 'Q2', title = 'Percent Cover for Caribbean Reefs (Median)')

allSummaries <- getSummaryStats(subset(df, (df$YEAR > 2010) & df$YEAR < 2015), substrateNames)
substrateBarplot(allSummaries, 'Q2', title = 'Percent Cover for Caribbean Reefs (Median)')

allSummaries <- getSummaryStats(subset(df, (df$YEAR > 2015) & df$YEAR < 2020), substrateNames)
substrateBarplot(allSummaries, 'Q2', title = 'Percent Cover for Caribbean Reefs (Median)')

###############
###############
# Clustering epxloration
###############
###############
# create train test split
train_ind <- sample(seq_len(nrow(df)), size = floor(0.75 * nrow(df)))
train <- df[train_ind, ]
test <- df[-train_ind, ]

# get new and old subset for testing
substrateNames <- subsetCols('substrate')
substrateNames <- c('HC','NI','RC','SD','RB','SC','SP')
#substrateNames <- c('HC','NI')
new <- test[as.Date(test$DATE) > as.Date('2010-06-01'),substrateNames]
old <- test[as.Date(test$DATE) <= as.Date('2010-06-01'),substrateNames]
x <- na.omit(train[,substrateNames])

# limit to correct colsj
#x <- na.omit(df[,substrateNames[c(1,6)]])
train <- na.omit(train[,substrateNames])
test <- na.omit(test[,substrateNames])
allSummaries <- getSummaryStats(df, substrateNames)
rownames(allSummaries) <- 1:length(allSummaries[,1])

# K Means
k <- kmeans(x, centers=4)
k$centers

# get percent changes from mean
percentChanges <- k$centers
for (i in 1:dim(percentChanges)[1]) {
  for (j in 1:length(substrateNames)) {
    meanIndex <- match(cols[j],allSummaries$NAME)
    percentChanges[i,j] <- 100*(percentChanges[i,j]-allSummaries$MEAN[meanIndex])/allSummaries$MEAN[meanIndex]
  }
}
print(percentChanges)

#' For both the Caribbean and non-Caribben datasets, there appear to be roughly 3 types of reefs (nore order changes on run)
#' 1. Healthy: high HC, SC, SP (and also SD and RB)
#' 2. Rocky: high RC low everything else
#' 3. Degraded: high NI low everything else

table(k$cluster)

######## test #########

# predict
table(predict.kmeans(k, old)) # ~43% healthy
table(predict.kmeans(k, new)) # ~18% healthy

# plot values (only good for 2D)
#plot(x, k$cluster)
#points(k$centers, col = 2, pch = 8, cex = 2)

# elbow method
k.max <- 15
data <- x
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

################
################
# Final
################
################
# setup
substrateNames <- c('HC','NI','RK','RC','RB','SD','SI','SC','SP')
set.seed(123)

# median barchart
allSummaries <- getSummaryStats(df, substrateNames)
substrateBarplot(allSummaries, 'Q2', title = 'Median Percent Cover for Caribbean Reefs by Substrate Type')

####### clustering ##########
# HC + NI, k = 2, all data (sanity check)
fit <- kmeans(na.omit(df[,c('HC','NI')]), centers=2, nstart = 25)
fviz_cluster(fit, data = na.omit(df[,c('HC','NI')]), main = 'KMeans Cluster for Caribbean Reef Data', xlab = 'Hard Coral',
             choose.vars = c('HC', 'NI'), 
             ylab = 'Nutrient Indicator Algae', geom = 'point', show.clust.cent = T) + theme_minimal() 

# RC, HC, NI, SD (4 biggest), k = 3, all data
out <- runKMeans(df = df, cols = c('HC','NI','RC','SD'), k =3, doTrainTest = F)
c1 <- substrateBarplot(out[1,]/100, 'C', title = '', isKMeans = T)
c2 <- substrateBarplot(out[2,]/100, 'C', title = '', isKMeans = T)
c3 <- substrateBarplot(out[3,]/100, 'C', title = '', isKMeans = T)
ggarrange(c3, c2,c1, 
          labels = c("1: Sand-Dominated Reefs", "2: Rock-Domindated Reefs", "3: Algae-Dominated Reefs"),
          ncol = 1, nrow = 3, hjust = -0.2)

# RC, HC, NI, SD, RB, SC, SP (4 biggest), k = 3, all data
out <- runKMeans(df = df, cols = c('HC','NI','RC','SD', 'RB', 'SC', 'SP'), k = 3, doTrainTest = F)
c1 <- substrateBarplot(out[1,]/100, 'C', title = '', isKMeans = T)
c2 <- substrateBarplot(out[2,]/100, 'C', title = '', isKMeans = T)
c3 <- substrateBarplot(out[3,]/100, 'C', title = '', isKMeans = T)
ggarrange(c3, c2,c1, 
          labels = c("1: Sand/Rubble-Dominated Reefs", "2: Rock-Domindated Reefs", "3: Algae-Dominated Reefs"),
          ncol = 1, nrow = 3, hjust = -0.2)

# 7 biggest substrate vals (proportionally), k = 4, all data
out <- runKMeans(df = df, cols = c('HC','NI','RC','SD', 'RB', 'SC', 'SP'), k = 4, doTrainTest = F)
c1 <- substrateBarplot(out[1,]/100, 'C', title = '', isKMeans = T)
c2 <- substrateBarplot(out[2,]/100, 'C', title = '', isKMeans = T)
c3 <- substrateBarplot(out[3,]/100, 'C', title = '', isKMeans = T)
c4 <- substrateBarplot(out[4,]/100, 'C', title = '', isKMeans = T)
ggarrange(c1, c4, c2,c3, 
          labels = c("1: Sand-Dominated Reefs", "2: Rock-Domindated Reefs", "3: Algae-Dominated Reefs",
                     "4: Coral-Dominated Reefs"),
          ncol = 2, nrow = 2, hjust = -0.2, font.label = list(size = 18))


# RC, HC, NI, SD, RB, SC, SP (4 biggest), k = 4, train test split (unused - inconsistent results)
#out <- runKMeans(df = df, cols = c('HC','NI','RC','SD', 'RB', 'SC', 'SP'), k = 4, doTrainTest = T, 
#                 splitDate = as.Date('2010-01-01'))
#substrateBarplot(out[1,]/100, 'C', title = 'Cluter 1: Sandy Reefs', isKMeans = T)
#substrateBarplot(out[2,]/100, 'C', title = 'Cluster 2: Nutrient Indicator Reefs', isKMeans = T)
#substrateBarplot(out[3,]/100, 'C', title = 'Cluster 3: Rocky Reefs', isKMeans = T)
#substrateBarplot(out[4,]/100, 'C', title = 'Cluster 3: Rocky Reefs', isKMeans = T)
