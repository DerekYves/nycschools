# Description: This script downloads tract-level US Census data for New York City and State
# Set these parameters only loaded when calling the script from outside the dio report process.


# User-defined Census parameters

# Refresh the Census Data?
refresh <- FALSE
acs_end_year  <- 2016 # Current ACS estimates
acs_cmpr_year <- 2010 # For retrieving decennial population estiamtes

report_statabbr <- "NY"

prjdir <- '~/projects/nycschools/data-raw/'

# Load libraries
library(rgdal)    # for readOGR and others
library(sp)       # for spatial objects
library(rgeos)    # For gintersection
library(leaflet)  # for interactive maps (NOT leafletR here)
library(ggplot2)  # for plotting
library(tigris)   # for tabular/shapefile merges
library(acs)      # for downloading census data
library(dplyr)    # run last due to ACS plyr dependency
library(maptools) # for fortify
library(scales)   # Pretty breaks in GGPlot

# install.packages(c("rgeos", "leaflet", "tigris", "acs", "maptools"))
# install.packages("tigris", type="source")

# library(stringr)  # to pad the FIPS codes

# Non-interactive report pamaters
key <- "YOUR KEY HERE"
geoparm <- geo.make(state=c(report_statabbr), county="*", tract="*")
api.key.install(key=key)

setwd(prjdir)

# Load the fips code for the state(s) to be mapped
fipmat <- read.csv(file="census/StateFIPSicsprAB.csv")
fips   <- sprintf("%02d", fipmat$FIPS[fipmat$AB %in% report_statabbr])


###### Census ACS Variables, 2016, 5 Year ######
# Varible Descriptions: http://api.census.gov/data/2014/acs5/variables.html

# See also:  http://api.census.gov/data/2009/acs5/variables.html


## (1) Population for the Geographic Unit ##
# B01003_001E	Total	B01003. Total Population [this is the population estimate]
varlist_pop <- "B01003_001"

## (2) Sample Size for the Geographic Unit ##

# *** NOT AVAILABLE FOR 2016 ACS***
# B00001_001E	Total	B00001. Unweighted Sample Count of the Population
# varlist_sample_size <- "B00001_001"


## (3) Income Variables ##
# Table: B19013. Median Household Income in the Past 12 Months (in 2014 Inflation-Adjusted Dollars) [this is simply the first column of the previous data pull. not used]
varlist_med_inc <- "B19019_001"

# Table: B17001. Poverty Status in the past 12 Months by Sex by Age
# -Variables-
# B17001_001E:	Total
# B17001_002E	Income in the past 12 months below poverty level
varlist_pov <- sprintf("B17001_%03d", 1:2)

# Table: B19001. Household Income in the Past 12 Months (in 2014 Inflation-Adjusted Dollars)
# B19001_001E Total:
# B19001_002E Less than $10,000
# B19001_003E $10,000 to $14,999
# B19001_004E $15,000 to $19,999
# B19001_005E $20,000 to $24,999
# B19001_006E $25,000 to $29,999
# B19001_007E $30,000 to $34,999
# B19001_008E $35,000 to $39,999
# B19001_009E $40,000 to $44,999
# B19001_010E $45,000 to $49,999
# B19001_011E $50,000 to $59,999
# B19001_012E $60,000 to $74,999
# B19001_013E $75,000 to $99,999
# B19001_014E $100,000 to $124,999
# B19001_015E $125,000 to $149,999
# B19001_016E $150,000 to $199,999
# B19001_017E $200,000 or more
varlist_hinc <- sprintf("B19001_%03d", 1:17)
var_inc_100k_up <- paste0(sprintf("B19001_%03d", c(14:17)), "$", collapse = "|")
var_inc_150k_up <- paste0(sprintf("B19001_%03d", c(16:17)), "$", collapse = "|")

## (4) Race ##
# Table: B02001. Race
# -Variables-
# B02001_001E Total:
# B02001_002E White alone
# B02001_003E Black or African American alone
# B02001_004E American Indian and Alaska Native alone
# B02001_005E Asian alone
# B02001_006E Native Hawaiian and Other Pacific Islander alone
# B02001_007E Some other race alone
# B02001_008E Two or more races:
# B02001_009E Two or more races:!!Two races including Some other race
# B02001_010E Two or more races:!!Two races excluding Some other race, and three or more races
varlist_race <- sprintf("B02001_%03d", 1:10)

var_non_wht <- paste0(sprintf("B02001_%03d", c(3:8)), "$", collapse = "|")

# Table: B03001. Hispanic or Latino by Specific Origin
# B03001_001E Total:
# B03001_002E Not Hispanic or Latino
# B03001_003E Hispanic or Latino
varlist_latin <- sprintf("B03001_%03d", 1:3)

## (5) Health Insurance Coverage ##
# Table: B27019 HEALTH INSURANCE COVERAGE STATUS AND TYPE BY AGE BY EDUCATIONAL ATTAINMENT
# -Variables-
# B27019_001E	Total:	B27019.
# B27019_002E	25 to 64 years:	B27019. [use this total]
# B27019_003E	25 to 64 years:!!Less than high school graduate:	B27019.
# B27019_004E	25 to 64 years:!!Less than high school graduate:!!With health insurance coverage	B27019.
# B27019_005E	25 to 64 years:!!Less than high school graduate:!!With health insurance coverage!!With private health insurance	B27019.
# B27019_006E	25 to 64 years:!!Less than high school graduate:!!With health insurance coverage!!With public coverage	B27019.
# B27019_007E	25 to 64 years:!!Less than high school graduate:!!No health insurance coverage	B27019.
# B27019_008E	25 to 64 years:!!High school graduate (includes equivalency):	B27019.
# B27019_009E	25 to 64 years:!!High school graduate (includes equivalency):!!With health insurance coverage	B27019.
# B27019_010E	25 to 64 years:!!High school graduate (includes equivalency):!!With health insurance coverage!!With private health insurance	B27019.
# B27019_011E	25 to 64 years:!!High school graduate (includes equivalency):!!With health insurance coverage!!With public coverage	B27019.
# B27019_012E	25 to 64 years:!!High school graduate (includes equivalency):!!No health insurance coverage	B27019.
# B27019_013E	25 to 64 years:!!Some college or associate's degree:	B27019.
# B27019_014E	25 to 64 years:!!Some college or associate's degree:!!With health insurance coverage	B27019.
# B27019_015E	25 to 64 years:!!Some college or associate's degree:!!With health insurance coverage!!With private health insurance	B27019.
# B27019_016E	25 to 64 years:!!Some college or associate's degree:!!With health insurance coverage!!With public coverage	B27019.
# B27019_017E	25 to 64 years:!!Some college or associate's degree:!!No health insurance coverage	B27019.
# B27019_018E	25 to 64 years:!!Bachelor's degree or higher:	B27019.
# B27019_019E	25 to 64 years:!!Bachelor's degree or higher:!!With health insurance coverage	B27019.
# B27019_020E	25 to 64 years:!!Bachelor's degree or higher:!!With health insurance coverage!!With private health insurance	B27019.
# B27019_021E	25 to 64 years:!!Bachelor's degree or higher:!!With health insurance coverage!!With public coverage	B27019.
# B27019_022E	25 to 64 years:!!Bachelor's degree or higher:!!No health insurance coverage	B27019.
varlist_hic <- sprintf("B27019_%03d", 1:22)

var_any_hc_u65 <- paste(sprintf("B27019_%03d", seq(4, 19, 5)), collapse = "|")
var_prv_hc_u65 <- paste(sprintf("B27019_%03d", seq(5, 20, 5)), collapse = "|")
var_pub_hc_u65 <- paste(sprintf("B27019_%03d", seq(6, 21, 5)), collapse = "|")
var_no_hc_u65  <- paste(sprintf("B27019_%03d", seq(7, 23, 5)), collapse = "|")

# Table: C27012. Health Insurance Coverage Status and Type by Work Experience by Age
# -Variables-
# C27012_003E	Worked full-time, year-round:!!18 to 64 years:
# C27012_004E	Worked full-time, year-round:!!18 to 64 years:!!With health insurance coverage
varlist_hic_emp <- sprintf("C27012_%03d", 3:4)

## (6) Education ##
# Table: B15003, EDUCATIONAL ATTAINMENT FOR THE POPULATION 25 YEARS AND OVER
# -Variables-
# B15003_001 Total
# B15003_002 No schooling completed
# B15003_003 Nursery school
# B15003_004 Kindergarten
# B15003_005 1st grade
# B15003_006 2nd grade
# B15003_007 3rd grade
# B15003_008 4th grade
# B15003_009 5th grade
# B15003_010 6th grade
# B15003_011 7th grade
# B15003_012 8th grade
# B15003_013 9th grade
# B15003_014 10th grade
# B15003_015 11th grade
# B15003_016 12th grade, no diploma
# B15003_017 Regular high school diploma
# B15003_018 GED or alternative credential
# B15003_019 Some college, less than 1 year
# B15003_020 Some college, 1 or more years, no degree
# B15003_021 Associate's degree
# B15003_022 Bachelor's degree
# B15003_023 Master's degree
# B15003_024 Professional school degree
# B15003_025 Doctorate degree
varlist_ed <- sprintf("B15003_%03d", c(1:25))

var_ed_no_hsch <- paste0(sprintf("B15003_%03d", c(02:16)), "$", collapse = "|")
var_ed_col_deg <- paste0(sprintf("B15003_%03d", c(22:25)), "$", collapse = "|")

## Functions to build variables
cen_sum <- function(df, var_group) {
    sum <- apply(df[ , grep(var_group, names(df))], 1, sum)
    return(sum)
}

cen_pct <- function(df, numerator, denominator, med = FALSE){
    if(med==T){
        percent  <- round(100 * (df[[numerator]] / median(df[[denominator]], na.rm=T)), digits = 1)
    } else {
        percent  <- round(100 * (df[[numerator]] / df[[denominator]]), digits = 1)
    }
    return(percent)
}

cen_ppct <- function(var){
    ppt <- paste0(var, "%")
    return(ppt)
}

# Build a vector of variable names
allvar <- as.vector(unlist(mget(ls(pattern = "varlist_"))))

## BUILD DATA FRAME

if(refresh) {
    # Download census data using the ACS/geoparm parameters set above
    acs_raw <- acs.fetch(endyear = acs_end_year, span = 5, geo = geoparm, var = allvar)

    # Now download the comparator population data from the 2010 decennial Census
    acs_prev_raw <- acs.fetch(endyear = acs_cmpr_year, span = 0, geo = geoparm, var = "P0010001", dataset="sf1")

    ## Save an image of Census Data for reuse
    save(acs_raw, acs_prev_raw, file=paste0(prjdir, "census/", report_statabbr, "_ACS5_", acs_end_year, ".rdata"))


} else {
    load(file= paste0(prjdir, "census/", report_statabbr, "_ACS5_", acs_end_year, ".rdata") )
}


# View the data frames
# if(interactive()) View(acs_raw@estimate)
# if(interactive() ) View(acs_prev_raw@estimate)


# Function to convert ACS class object to a data.frame to simplify merging
acs_trans <- function(df) {
    x <- data.frame(
        GEOID = paste0(
            sprintf('%02d', df@geography$state),
            sprintf('%03d', df@geography$county),
            sprintf('%06d', as.numeric(df@geography$tract))),
        state=sprintf('%02d', df@geography$state),
        county=sprintf('%03d', df@geography$county),
        tract=sprintf('%06d', as.numeric(df@geography$tract)),
        description=rownames(df@estimate),
        df@estimate,
        stringsAsFactors = FALSE,
        row.names = NULL)
    return(x)
}

acs_df      <- acs_trans(acs_raw)
acs_df_prev <- acs_trans(acs_prev_raw)

# Rename the previous year's population variable and join with current

popframe <- select_(acs_df_prev, "GEOID", decennial_10_pop="P0010001")
acs_df   <- left_join(acs_df, popframe, by="GEOID", all.x=T)

# Tracts with missing 2010 population data
table(is.na(acs_df$decennial_10_pop))

### Perform Variable Transformations ###

## Population ##
# Population Change (currently '09-'14 ACS estiamte comparisons)
acs_df$pct_pop_change <- with(acs_df,
                              ifelse(!is.na(B01003_001) & !is.na(decennial_10_pop) & decennial_10_pop != 0,
                                     round(100 * ((B01003_001 - decennial_10_pop) / decennial_10_pop), digits=1), NA))

range(acs_df$pct_pop_change, na.rm=T)
acs_df$ppct_pop_change <- cen_ppct(acs_df$pct_pop_change)

# Recode numeric percent population change to +/1 50%
acs_df$pct_pop_change[acs_df$pct_pop_change >= 15] <- 15
acs_df$pct_pop_change[acs_df$pct_pop_change <= -15] <- -15

## Income ##

# Poverty rate (B17001_001)
acs_df$pct_pov <- cen_pct(acs_df, "B17001_002", "B17001_001")
acs_df$ppct_pov <- cen_ppct(acs_df$pct_pov)
acs_df$pct_pov <- cen_pct(acs_df, "B17001_002", "B17001_001")
acs_df$pct_pov[acs_df$pct_pov >= 50] <- 50

summary(acs_df$pct_pov)

# Income 100k and up (B19001_001E)
acs_df$sum_inc_100k_up  <- cen_sum(acs_df, var_inc_100k_up)
acs_df$pct_inc_100k_up  <- cen_pct(acs_df, "sum_inc_100k_up", "B19001_001")
acs_df$ppct_inc_100k_up <- cen_ppct(acs_df$pct_inc_100k_up)
acs_df$pct_inc_100k_up[acs_df$pct_inc_100k_up >= 50] <- 50

# Income 150k and up (B19001_001E)
acs_df$sum_inc_150k_up  <- cen_sum(acs_df, var_inc_150k_up)
acs_df$pct_inc_150k_up  <- cen_pct(acs_df, "sum_inc_150k_up", "B19001_001")
acs_df$ppct_inc_150k_up <- cen_ppct(acs_df$pct_inc_150k_up)

# Ratio of median income to the catchment area
# Missing values are coded as '-666666666'

table( acs_df[['B19019_001']][acs_df$B19019_001 <1 ])
acs_df[['B19019_001']][acs_df$B19019_001 <1 ] <- NA
summary(acs_df[['B19019_001']])
acs_df$est_med_hinc  <- acs_df[['B19019_001']]
acs_df$ratio_med_hinc  <- cen_pct(acs_df, "B19019_001", "B19019_001", med=T)
acs_df$pratio_med_hinc <- cen_ppct(acs_df$ratio_med_hinc)

## Race/Ethnicity ##

# Percentages categories ("alone" or two or more)(B02001_001)
acs_df$pct_white  <- cen_pct(acs_df, "B02001_002", "B02001_001")
acs_df$ppct_white <- cen_ppct(acs_df$pct_white)

acs_df$pct_black  <- cen_pct(acs_df, "B02001_003", "B02001_001")
acs_df$ppct_black <- cen_ppct(acs_df$pct_black)

acs_df$pct_amind  <- cen_pct(acs_df, "B02001_004", "B02001_001")
acs_df$ppct_amind <- cen_ppct(acs_df$pct_amind)

acs_df$pct_asian  <- cen_pct(acs_df, "B02001_005", "B02001_001")
acs_df$ppct_asian <- cen_ppct(acs_df$pct_asian)

acs_df$pct_nhpi   <- cen_pct(acs_df, "B02001_006", "B02001_001")
acs_df$ppct_nhpi  <- cen_ppct(acs_df$pct_nhpi)

acs_df$pct_oth    <- cen_pct(acs_df, "B02001_007", "B02001_001")
acs_df$ppct_oth   <- cen_ppct(acs_df$pct_oth)

acs_df$pct_mlt_r  <- cen_pct(acs_df, "B02001_008", "B02001_001")
acs_df$ppct_mlt_r <- cen_ppct(acs_df$pct_mlt_r)

# Percent of non-white residents and ratio to the diocesan catchment area
acs_df$sum_nonwht  <- cen_sum(acs_df, var_non_wht)
acs_df$pct_nonwht  <- cen_pct(acs_df, "sum_nonwht", "B02001_001")
acs_df$ppct_nonwht <- cen_ppct(acs_df$pct_nonwht)

# Latino percentage categories (B03001_001)
acs_df$pct_latin  <- cen_pct(acs_df, "B03001_003", "B03001_001")
acs_df$ppct_latin <- cen_ppct(acs_df$pct_latin)

## Health Insurance ##

# Coverage types for under 65 population (B27019_002)
acs_df$sum_any_hc_u65  <- cen_sum(acs_df, var_any_hc_u65)
acs_df$pct_any_hc_u65  <- cen_pct(acs_df, "sum_any_hc_u65", "B27019_002")
acs_df$ppct_any_hc_u65 <- cen_ppct(acs_df$pct_any_hc_u65)

acs_df$sum_prv_hc_u65  <- cen_sum(acs_df, var_prv_hc_u65)
acs_df$pct_prv_hc_u65  <- cen_pct(acs_df, "sum_prv_hc_u65", "B27019_002")
acs_df$ppct_prv_hc_u65 <- cen_ppct(acs_df$pct_prv_hc_u65)

acs_df$sum_pub_hc_u65  <- cen_sum(acs_df, var_pub_hc_u65)
acs_df$pct_pub_hc_u65  <- cen_pct(acs_df, "sum_pub_hc_u65", "B27019_002")
acs_df$ppct_pub_hc_u65 <- cen_ppct(acs_df$pct_pub_hc_u65)

acs_df$sum_no_hc_u65   <- cen_sum(acs_df, var_no_hc_u65)
acs_df$pct_no_hc_u65   <- cen_pct(acs_df, "sum_no_hc_u65", "B27019_002")
acs_df$ppct_no_hc_u65  <- cen_ppct(acs_df$pct_no_hc_u65)

## Education Variables ##

# College degree or higher (B15003_001)
acs_df$sum_unv_deg   <- cen_sum(acs_df, var_ed_col_deg)
acs_df$pct_unv_deg   <- cen_pct(acs_df, "sum_unv_deg", "B15003_001")
acs_df$ppct_unv_deg  <- cen_ppct(acs_df$pct_unv_deg)

# Less than high school
acs_df$sum_under_hs  <- cen_sum(acs_df, var_ed_no_hsch)
acs_df$pct_under_hs  <- cen_pct(acs_df, "sum_under_hs", "B15003_001")
acs_df$ppct_under_hs <- cen_ppct(acs_df$pct_under_hs)

# Recode low count/missing GEOIDs to NA
if(exists("report_empty_geoids")){
    acs_df[acs_df$GEOID %in% report_empty_geoids, -1] <- NA
}

## Clean the dataframe
acs_df[acs_df=="NaN"]  <- NA
acs_df[acs_df=="NaN%"] <- NA
acs_df[acs_df=="NA%"]  <- NA
acs_df[acs_df=="Inf"]  <- NA
acs_df[acs_df=="Inf%"] <- NA

names(acs_df)

## Reorder the dataframe
acs <- cbind(GEOID=acs_df$GEOID,
             state=acs_df$state,
             county=acs_df$county,
             tract=acs_df$tract,
             description=acs_df$description,
             acs_df[ , grep("pct_|sum_|^est_|ratio_", names(acs_df))],
             acs_df[ , grep("^B", names(acs_df))], stringsAsFactors=FALSE)




## Data Quality
# DQ:  Pecentages should have a max of 100, with no negative numbers
summary(acs[ , grep("^pct_", names(acs))])
# DQ:  Ratios should have no negative numbers
summary(acs[ , grep("^ratio_", names(acs))])
# Income should be positive
summary(acs[ , grep("^est_", names(acs))])

## Save reduced frame
acsr <- select(acs, -starts_with("B"), -starts_with("ppct"))
saveRDS(acsr, paste0(prjdir, "census/", report_statabbr, "_ACS5_df_", acs_end_year, ".Rds"))
