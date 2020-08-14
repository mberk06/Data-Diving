#' Author: Michael Berk
#' Date: August 2020
#' Descrtiption: 
#'     This script was used to perform a visual analysis on TS data in the Caribbean. The data are
#'     supplied by Reef Check, a non-profit that trains amatuer daivers around the world to collect
#'     data. Data limited to the Caribbean includes 1576 unique dives from 1997-05-24 to 2019-08-24.
#'     
#'     The file is roughly broken down into three sections:
#'     1. Helpers (which are called by subsequent code)
#'     2. Exploration. This code is messy and was used to inform analysis direction.  
#'     3. TS Plots. These were the plots that were used in the post. However, some plots are commented
#'        out or marked as "unused" because their conclusions were deemed repetetive or unrelated to
#'        the post. 
#'     
#'     The steps are as follows:
#'     1. Aggregate data annually using a mean agg func. Annual timeframes were selected because 
#'        they account for seasonality. In future iterations, quarterly or bi-annual timeframes 
#'        may produce interesting results.
#'     2. The data were plotted and a weighted linear regression was fit. Weights were chosen to handle 
#'        an uneven dive distribution (the tails of the TS exhibited less dives, especially for the 
#'        last 5 years); each wieght corresponds to the number of dives that year.
#'     
#'     Notes:
#'     - Not all organisms experienced intuitive trends. For instance, in 2018 and 2019, many fish 
#'       species exhibited extremely high levels. Anomolies observed here will be tackled in future 
#'       analyses.
#'     - Rock and diadema showed strong positive trends.
#'     - Further areas of analysis are listed outside this documentation.
#'     
#'     Feel free to reach out through github or at https://michaeldberk.com/contact.

pacman::p_load(ggplot2) 
pacman::p_load(gghighlight)
pacman::p_load(lubridate) # truncate dates
pacman::p_load(stringr) # str_replace
pacman::p_load(dplyr) # aggregate counts


# set wd based on user
if (dir.exists("/Users/michaelberk")) {
  setwd('~/Documents/Ocean Analysis/')
} else {
  setwd('~ISA PATH') # TODO CHANGE TO CORRECT PATH
}

# setup data
source(paste0(getwd(), '/Scripts/setupData.R'))

#######################
#######################
# Helpers
#######################
#######################
aggDF <- function(df, y, timeFrame) {
  #' aggreagte df for timeframe
  #' 
  #' @param df data.frame to be aggregated (data.frame)
  #' @param y y variable col name to be modeled (character)
  #' @param timeframe determine the interval to agg data (character)
  #' @return aggregated df for given timeframe (data.frame)
  
  # parse timeframe
  if (timeFrame == 'DATE') {
    df$DATE <- as.Date(df$DATE)
  } else if (timeFrame == 'WEEK') {
    df$WEEK <- round_date(as.Date(df$DATE), 'week')
  } else if (timeFrame == 'MONTH') {
    df$MONTH <- round_date(as.Date(df$DATE), 'month') 
  } else if (timeFrame == 'QUARTER') {
    df$QUARTER <- round_date(as.Date(df$DATE), '3 months')
  } else if (timeFrame == 'BIANNUAL') {
    df$BIANNUAL <- round_date(as.Date(df$DATE), '6 months')
  } else if (timeFrame == 'YEAR') {
  } else {
    print('Incorrect timeframe')
    return()
  }
  
  # used for modeling
  #cols <- c(subsetCols('organism'), timeFrame, subsetCols('substrate'), subsetCols('anthro'))
  #dfToPivot <- df[,cols]
  
  # create count df to pivot
  dfToPivot <- df %>% count(get(timeFrame), get(y))
  names(dfToPivot) <- c(timeFrame, y, 'N')
  
  # aggregate df on year
  aggDFOrganism <- aggregate(get(y) ~ get(timeFrame), data = dfToPivot, FUN = mean)
  aggDFCount <- aggregate(N ~ get(timeFrame), data = dfToPivot, FUN = sum)
  aggDF <- merge(aggDFOrganism, aggDFCount, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE)
  
  # rename cols
  names(aggDF) <- c('DATE', y, 'N')
  
  # convert Date (if necessary) and sort
  if (timeFrame != 'YEAR') {
    aggDF$DATE <- as.Date(aggDF$DATE)
  }
  aggDF <- aggDF[order(aggDF$DATE),]
  
  # return
  return(aggDF)
}

subsetCols <- function(subsetType) {
  #' return column names for specific subset
  #' 
  #'  @param subsetType key of list that is used to return col name subset (character)
  #'  @return vector of column names (character)
  
  subsets <- list()
  subsets[['substrate']] <- c("FS","HC","NI","OT","RB","RC","RK","SC","SD","SI","SP")
  subsets[['anthro']] <- c("DYNAMITE_FISHING","POISON_FISHING","AQUARIUM_FISH_COLLECTION","HARVEST_OF_INVERTS_FOR_FOOD",
                           "HARVEST_OF_INVERTS_FOR_CURIO","TOURIST_DIVING_SNORKELING","SEWAGE_POLLUTION",
                           "INDUSTRIAL_POLLUTION",
                           'SILTATION','RECREATIONAL_FISHING','INVERTEBRATE_SHELL_COLLECTION', 'ANCHORING', 
                           'DIVING',"TRASH_FISH_NETS","TRASH_GENERAL","CORAL_DAMAGE_ANCHOR","CORAL_DAMAGE_DYNAMITE",
                           "CORAL_DAMAGE_OTHER","SPEARFISHING") 
  subsets[['organismAll']] <- c('ARABIAN_BUTTERFLYFISH','ASPERGILLOSIS','BANDED_CORAL_SHRIMP','BARRACUDA',
                                'BARRAMUNDI_COD','BLACK_BAND','BLACK_BAND','BLACK_SPOTTED_GRUNT',
                                'BLACK_URCHIN','BLEACHING_PER_OF_COLONY','BLEACHING_PER_OF_POPULATION',
                                'BLUE_SEA_STAR','BLUELINE_SNAPPER','BROOMTAIL_WRASSE',
                                'BUMPHEAD_PARROT','BUTTERFLYFISH','COTS',
                                'COWRIES','DARK_BUTTERFLYFISH','DARK_BUTTERFLYFISH','DIADEMA',
                                'EDIBLE_SEA_CUCUMBER','FLAMINGO_TONGUE','GIANT_CLAM_10-20_CM',
                                'GIANT_CLAM_20-30_CM','GIANT_CLAM_30-40_CM','GIANT_CLAM_40-50_CM',
                                'GIANT_CLAM_<10_CM','GIANT_CLAM_>50_CM','GIANT_CLAM_TOTAL',
                                'GIANT_HAWKFISH','GOATFISH','GORGONIAN','GREY_GRUNT','GROUPER_30-40_CM',
                                'GROUPER_40-50_CM','GROUPER_50-60_CM','GROUPER_>60_CM','GROUPER_TOTAL',
                                'GRUNTS','HAEMULIDAE','HELMET_CONCH','HUMPHEAD_WRASSE',
                                'JACKS','KING_ANGELFISH','LIONFISH','LOBSTER','LONGFIN_BANNERFISH',
                                'MANTAS','MEXICAN_HOGFISH','MORAY_EEL','NASSAU_GROUPER_30-40_CM',
                                'NASSAU_GROUPER_40-50_CM','NASSAU_GROUPER_50-60_CM',
                                'NASSAU_GROUPER_>60_CM','NASSAU_GROUPER_TOTAL','NASSAU_GROUPER_30-40_CM',
                                'NASSAU_GROUPER_40-50_CM','NASSAU_GROUPER_50-60_CM',
                                'NASSAU_GROUPER_GREATER_60_CM','NASSAU_GROUPER_TOTAL','ORANGE_SPINE_UNICORNFISH',
                                'ORANGE_SPOTTED_GROUPER_30-40_CM','ORANGE_SPOTTED_GROUPER_40-50_CM',
                                'ORANGE_SPOTTED_GROUPER_50-60_CM','ORANGE_SPOTTED_GROUPER_>60_CM',
                                'ORANGE_SPOTTED_GROUPER_TOTAL','PARROTFISH',
                                'PEACOCK_GROUPER_30-40_CM','PEACOCK_GROUPER_40-50_CM',
                                'PEACOCK_GROUPER_50-60_CM','PEACOCK_GROUPER_GREATER_60_CM',
                                'PEACOCK_GROUPER_TOTAL','PENCIL_URCHIN','QUEEN_CONCH','SEA_FAN','SHARKS',
                                'SHORT_SPINE_URCHIN','SLATE_PENCIL_URCHIN','SNAPPER',
                                'SPIDER_CRAB','SPOTTED_GRUNT','TRIPNEUSTES','TRITON','TROCHUS','TURTLES',
                                'WHITE_BAND','WHITE_PLAGUE','WHITE_BAND','YELLOW_GOATFISH','YELLOW_TANG',
                                'YELLOWBAR_ANGELFISH','YELLOWTAIL_TANG')
  subsets[['organism']] <- c('SNAPPER','TRIPNEUSTES','TRITON','PENCIL_URCHIN','PARROTFISH','MORAY_EEL','LOBSTER','GROUPER_TOTAL',
                             'DIADEMA','BUTTERFLYFISH','BANDED_CORAL_SHRIMP','BLEACHING_PERCENT_OF_COLONY','HAEMULIDAE',
                             'BLEACHING_PERCENT_OF_POPULATION') 
  subsets[['percentile']] <- unlist(lapply(names(df), function(x) x[grepl("PERCENTILE", x)]))
  
  # create factors for ordering
  subsets[['f1']] <- c('ANCHORING','DIVING','RECREATIONAL_FISHING','INVERTEBRATE_SHELL_COLLECTION','SPEARFISHING')# yes/no factors
  subsets[['f2']] <- c('SILTATION','') # ALWAYS, OFTEN, OCCASIONALLY, NEVER
  subsets[['f3']] <- c('HARVEST_OF_INVERTS_FOR_CURIO','TOURIST_DIVING_SNORKELING', # HIGH, MODERATE, LOW, NONE
                       'SEWAGE_POLLUTION','INDUSTRIAL_POLLUTION','DYNAMITE_FISHING','POISON_FISHING','AQUARIUM_FISH_COLLECTION',
                       'HARVEST_OF_INVERTS_FOR_FOOD')
  
  # return if in susbets
  if (subsetType %in% names(subsets)) {
    subsets[[subsetType]] <- toupper(gsub("/","\\.",gsub("\\?","\\.",gsub(" ","\\.",subsets[[subsetType]]))))
    subsets[[subsetType]] <- ifelse(substr(subsets[[subsetType]],nchar(subsets[[subsetType]]),nchar(subsets[[subsetType]]))==".",substr(subsets[[subsetType]],1,nchar(subsets[[subsetType]])-1),subsets[[subsetType]])
    subsets[[subsetType]] <- str_replace(subsets[[subsetType]], "[-]", ".")
    subsets[[subsetType]] <- str_replace(subsets[[subsetType]], "[>]", "")
    subsets[[subsetType]] <- str_replace(subsets[[subsetType]], "[<]", "")
    return(subsets[[subsetType]]) 
  } else {
    print("Incorrect subset key.")
    return(NA)
  }
}

tsPlot <- function(y, title, xcoordOffset, ycoordDivisor, isCoverage) {
  #' Create ggplot TS plot with yearly aggregtion for a given column
  #' 
  #' @param y (character) col name to aggregate and plot
  #' @param title (character) title of col name to go in ggtitle
  #' @param xcoordOffset (numeric) positive or negative number to shift the slope label
  #' @param ycoordDivisor (numeric) value that determines how far above the plotted lm to label the slope
  #' @param isCoverage (bool) if is % coverage in title and axis label
 
  # pars isCoverge
  if (isCoverage) {
    ylab = '% Coverage'
    title2 = ' % Coverage for Caribbean Reefs'
  } else {
    ylab = 'Count Per Dive'
    title2 = ' Counts for Caribbean Reefs'
  }
  
  # create df
  adf = aggDF(df, y, timeFrame = 'YEAR')
  names(adf) <- c('DATE', 'Y', 'N')
  
  # convert to percent
  if (y %in% subsetCols('substrate')) {
    adf[,2] <- adf[,2]*100
  }
  
  # create lm
  lm <- lm(Y ~ DATE, data = adf, weights = N)
  yint <- lm$coefficients[[1]]
  slope <- lm$coefficients[[2]]
  xcoord <- min(adf$DATE) + (max(adf$DATE) - min(adf$DATE))/2 + xcoordOffset
  predicted_df <- data.frame(y_hat = predict(lm, adf), DATE = adf$DATE)
  
  print('lm and fitted values')
  print(lm)
  print(predicted_df)
  
  # create plot 
  p <- ggplot(data = adf, aes(x = DATE, y = Y)) + 
    geom_line(color = blue, size = 1) +
    theme_minimal() + xlab('') + ylab(ylab) +
    ggtitle(paste0(title, title2)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = .8), axis.text.y = element_text(size=10, face='bold')) +
    scale_x_continuous(breaks = round(seq(min(adf$DATE), max(adf$DATE), by = 2),1)) +
    geom_line(color=red, data = predicted_df, aes(x=DATE, y=y_hat))
  
    
  # create ranges for slope label
  y_range <- layer_scales(p)$y$range$range
  x_range <- layer_scales(p)$x$range$range
  x_to_y <- (x_range[2] - x_range[1])/(y_range[2] - y_range[1])
  
  # create label
  p <- p + annotate('text', x = xcoord, y = xcoord * slope + yint*1 + mean(adf[,2])/ycoordDivisor,
             angle = atan(slope * x_to_y) * 180/pi, 
             label = paste0('m = ', round(slope, digits = 3)))
  p
}

highlightedRegion <- function(p, x1, x2) {
  #' Create a highlighed region on a given plot starting at x1 and ending at x2\
  #' 
  #' @param p (ggplot plot) plot to add to
  #' @param x1 (numeric or date) x coordinate to start highlight. Can be vector.
  #' @param x2 (numeric or date) x coordinate to end highlight. Can be vector.
  #' @return ggplot plot
  
  # get y range
  y_range <- layer_scales(p)$y$range$range
  
  # create df with values 
  rectDF <- data.frame(x1 = x1, x2 = x2)
  
  # create rect and add to plot
  p <- p + geom_rect(data=rectDF, inherit.aes=FALSE, 
            aes(xmin=x1, xmax=x2, ymin=y_range[1], ymax=y_range[2]), 
            color="transparent", fill=red, alpha=0.2)
}



#######################
#######################
# Exploration
#######################
#######################

plot(df$DATE, df$BUTTERFLYFISH)
lm(BUTTERFLYFISH ~ DATE, data = df)

for (o in subsetCols('organism')) {
  adf = aggDF(df, o, timeFrame = 'YEAR')
  plot(adf$DATE, adf[,2], type='l', main = o)
}

# create df
for (c in caribbeanCountries) {
  # subset to country and aggregate annually 
  if (c %in% df$COUNTRY) {
    temp <- subset(df, COUNTRY == c)
    adf <- aggDF(temp, 'GROUPER_TOTAL', timeFrame = 'YEAR')
    
    # get country change 2015 - 2018 (where there seems to be a big spike)
    change <- 100*(adf[length(adf$DATE) - 2,2] - adf[length(adf$DATE) - 5,2])/adf[length(adf$DATE) - 5,2]
    
    print('-----------------')
    print(paste0(adf[length(adf$DATE) - 5,2], paste0(' - ', adf[length(adf$DATE) - 2,2])))
    
    if (length(change) > 0) {
      print(paste0(paste0(c, ': '), change))
    }
  }
}
adf = aggDF(df, 'PARROTFISH', timeFrame = 'YEAR')
lm(PARROTFISH ~ DATE, data = adf, weights = N)

# why are grouper numbers so high
tempG <- subset(df, YEAR == 2018)
hist(tempG$GROUPER_TOTAL, breaks = 40)

#######################
#######################
# TS Plots
#######################
#######################
# bleaching timeframes (unused)
starts <- c(1997, 2005, 2015)
stops <- c(1998, 2006, 2017)
#p <- highlightedRegion(p, starts, stops)

# Hard coral
tsPlot(y = 'HC', title = 'Hard Coral', xcoordOffset = 2, ycoordDivisor = 40, isCoverage = T)

# Soft coral
tsPlot(y = 'SC', title = 'Soft Coral', xcoordOffset = -7.8, ycoordDivisor = 20, isCoverage = T)

# Rock (unused)
#tsPlot(y = 'RC', title = 'Rock', xcoordOffset = -2, ycoordDivisor = -35, isCoverage = T)

# butterfly fish
tsPlot(y = 'BUTTERFLYFISH', title = 'Butterfly Fish', xcoordOffset = 5.5, ycoordDivisor = 30, isCoverage = F)

# parrot fish
p <- tsPlot(y = 'PARROTFISH', title = 'Parrot Fish', xcoordOffset = 1, ycoordDivisor = 30, isCoverage = F)
p

# overlay HC (not used due to visual confusion)
#adf <- aggDF(df, 'HC', 'YEAR')
#p + geom_line(aes(x = adf[,1], y=adf[,2]*100)) 

# parrot fish
tsPlot(y = 'GROUPER_TOTAL', title = 'Grouper', xcoordOffset = 5.4, ycoordDivisor = 25, isCoverage = F)
