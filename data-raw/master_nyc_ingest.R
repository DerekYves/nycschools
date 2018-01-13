################################################################################
################################################################################
################################################################################

#### Master Script to Ingest NYC School Survey Data ####
rm(list=ls())

# !diagnostics off

# User Defined Parameters

na_strings <- c('N/A', '', 'Not Scored')
# prjdir <- '~/Dropbox/courses/soc345/data/nyc_school_surey/'

prjdir <- '~/projects/nycschools/data-raw/'
report_statabbr <- "NY"

years <- c(2013, 2016, 2017)

library(readr)
library(dplyr)
library(rgdal)    # for readOGR and others
library(proj4)



################################################################################
################################################################################
################################################################################

##### Load all files
setwd(prjdir)

#### 2013 Data
yr <- 1
cat("Starting", years[yr], "\n")

## NOTE: CSV files for this analysis were culled for the XLS files provided
# by the DOE. Due to heavy formatting and merged cells in the source, the CSV
# files read by this script were created manually. The original xls files are
# included in the raw-data folder and are named, e.g.: 2017PublicDataFileParent.xlsx

## 2013 Parent File (checked raw for q1c, 01M019, matches)
# Note: Data for this year has a composite score for each question ranging from 0-10

p <- read_csv(paste0(years[yr], '/data/parent_', years[yr], '_score.csv'), na=na_strings)

# View(p[, sapply(p, function(i) all(is.na(i)))])

p = p[, sapply(p, function(i) !all(is.na(i)))] # Drop entirely NA cols

# Should be empty
unique(select_if(select(p, -dbn, -locationname, -school_type), is.character))


## 2013 Student File (checked raw for q5a, 01M034/roosevelt, matches)
s <- read_csv(paste0(years[yr], '/data/student_', years[yr], '_score.csv'), na=na_strings)
s <- s[, sapply(s, function(i) !all(is.na(i)))]


## 2013 Total File (checked raw for 01M064, eng_tot_13, matches)
tot <- read_csv(paste0(years[yr], '/data/total_', years[yr], '_score.csv'), na=na_strings)
tot <- tot[, sapply(tot, function(i) !all(is.na(i)))]

# Should match
stopifnot(all(all.equal(nrow(s), nrow(p)), all.equal(nrow(s), nrow(tot))))


## Join the student/parent/toal files
test <- nrow(s)
out <- left_join(s, p)
out <- left_join(out, tot)
stopifnot(test==nrow(out))

# Write files to disk
out <- cbind(year=years[yr], out)
saveRDS(out, paste0(years[yr], '/data/merged_', years[yr], '_scores.Rds'))
write.csv(out, paste0(years[yr], '/data/merged_', years[yr], '_scores.Rds'), row.names=FALSE)

saveRDS(data.frame(dbn=unique(out$dbn)), paste0(years[yr], '/data/merged_', years[yr], '_dbns.Rds'))



################################################################################
################################################################################
################################################################################

#### 2016 Data
yr <- 2
cat("Starting", years[yr], "\n")

## 2016 Parent File: N=3,086 in raw source file (verified)
p <- read_csv(paste0(years[yr], '/data/parent_', years[yr], '_score.csv'), na=na_strings)
nrow(p)

# View(p[, sapply(p, function(i) all(is.na(i)))]) # only five empty cols, q5
p = p[, sapply(p, function(i) !all(is.na(i)))] # Drop entirely NA cols

p <- rename(p, p_rr_par=p_rr, locationname_par=locationname)

# Data Checks:  DBN 02M051
p$p_q1h_sa[p$dbn=='02M051'] # 100 (100 in nsource)
p$p_q2c_sd[p$dbn=='02M051'] # 10 (10 in source)
p$p_q5a_s[p$dbn=='02M051']  # 104 (104 in source)

# Should be empty
unique(select_if(select(p, -dbn, -locationname_par), is.character))


## 2016 Student File (N=1,137 in raw source file)
s <- read_csv(paste0(years[yr], '/data/student_', years[yr], '_score.csv'), na=na_strings)
nrow(s)

s <- s[, sapply(s, function(i) !all(is.na(i)))] # no change
s <- rename(s, s_rr_stu=s_rr, locationname_stu=locationname)

# Data Checks:  DBN 02M051
s$s_q1e_a[s$dbn=='02M418']     # 238 (238 in nsource)
s$s_q2d_most[s$dbn=='02M418']  # 312 (312in source)
s$s_q6e_none[ s$dbn=='02M418'] # 9 (9 in source)


## 2016 Total File: N=3,098 in raw source file
# (checked raw for 01M064, eng_tot_13, matches)
# Note: These contain a composite score, 0:1
tot <- read_csv(paste0(years[yr], '/data/total_', years[yr], '_score.csv'), na=na_strings)
nrow(tot)

tot <- tot[, sapply(tot, function(i) !all(is.na(i)))]
tot <- rename(tot, p_rr_tot=p_rr, s_rr_tot=s_rr, t_rr_tot=t_rr, locationname_tot=locationname)

tot$eff_leaders[tot$dbn=='05M030'] # .82 (82% in nsource)
tot$supp_env[tot$dbn=='05M030']    # .68 (68% in source)
tot$stu_sat[tot$dbn=='05M030']     #  NA (N/A in source)


## Join the parent/student/total files (inner on student)
test <- nrow(s)
out <- inner_join(s, p) # Reduce record count
nrow(s); (test<-nrow(out)) # lose no rows

out <- inner_join(tot, out) # lose no rows
stopifnot(test==nrow(out))

# These variables be similar or match exactly
# View(select(out, contains("rr")))
out$p_rr_par <- out$s_rr_stu <- NULL

# Write files to disk
out <- cbind(year=years[yr], out)
saveRDS(out, paste0(years[yr], '/data/merged_', years[yr], '_scores.Rds'))
write.csv(out, paste0(years[yr], '/data/merged_', years[yr], '_scores.csv'), row.names=FALSE)

saveRDS(data.frame(dbn=unique(out$dbn)), paste0(years[yr], '/data/merged_', years[yr], '_dbns.Rds'))

assign(paste0("surveys_", years[yr]), out)

# Stubs in this file
table(sub("^p_q\\d+._", "", names(p)))


################################################################################
################################################################################
################################################################################

#### 2017 Data
yr <- 3
cat("Starting", years[yr], "\n")

## 2017 Parent File: N=2,966 in raw source file (verified)
p <- read_csv(paste0(years[yr], '/data/parent_', years[yr], '_score.csv'), na=na_strings)
nrow(p)

# View(p[, sapply(p, function(i) all(is.na(i)))]) # only five empty cols, q5
p = p[, sapply(p, function(i) !all(is.na(i)))] # Drop entirely NA cols

p <- rename(p, p_rr_par=p_rr, locationname_par=locationname)

# Data Checks:  DBN 02M051
p$p_q1e_a[p$dbn== '02M051'] # 85  (85 in nsource)
p$p_q2b_sd[p$dbn=='02M051'] # 2   (2 in source)
p$p_q5a_vs[p$dbn=='02M051'] # 113 (113 in source)

# Should be empty
unique(select_if(select(p, -dbn, -locationname_par), is.character))

## 2017 Student File (N=1,121 in raw source file)
s <- read_csv(paste0(years[yr], '/data/student_', years[yr], '_score.csv'), na=na_strings)
nrow(s)

s <- s[, sapply(s, function(i) !all(is.na(i)))] # no change
s <- rename(s, s_rr_stu=s_rr, locationname_stu=locationname)

# Data Checks:  DBN 29Q248
s$s_q1h_a[s$dbn=='29Q248']     # 122 (122 in nsource)
s$s_q2d_most[s$dbn=='29Q248']  # 103 (103 in source)
s$s_q6e_none[ s$dbn=='29Q248'] # 50 (50 in source)


## 2017 Total File: N=3,254 in raw source file
# (checked raw for 01M064, eng_tot_13, matches)
# Note: These contain a composite score, 0:1
tot <- read_csv(paste0(years[yr], '/data/total_', years[yr], '_score.csv'), na=na_strings)
nrow(tot)

tot <- tot[, sapply(tot, function(i) !all(is.na(i)))]
tot <- rename(tot, p_rr_tot=p_rr, s_rr_tot=s_rr, t_rr_tot=t_rr, locationname_tot=locationname)

tot$eff_leaders[tot$dbn=='26Q191'] # 4.24 (4.24 in nsource)
tot$supp_env[tot$dbn=='26Q191']    # 4.14 (4.14in source)
tot$trust[ tot$dbn=='26Q191']      # 4.18 (4.18 in source)


## Join the parent/student/total files (inner on student)
test <- nrow(s)
out <- inner_join(s, p) # Reduce record count
nrow(s); (test<-nrow(out)) # lose 3 rows

out <- inner_join(tot, out) # lose no rows
stopifnot(test==nrow(out))

# These variables be similar or match exactly
# View(select(out, contains("rr")))
out$p_rr_par <- out$s_rr_stu <- NULL

# Write files to disk
out <- cbind(year=years[yr], out)
saveRDS(out, paste0(years[yr], '/data/merged_', years[yr], '_scores.Rds'))
write.csv(out, paste0(years[yr], '/data/merged_', years[yr], '_scores.csv'), row.names=FALSE)

assign(paste0("surveys_", years[yr]), out)
saveRDS(data.frame(dbn=unique(out$dbn)), paste0(years[yr], '/data/merged_', years[yr], '_dbns.Rds'))

# Stubs in this file
table(sub("^p_q\\d+._", "", names(p)))


################################################################################
################################################################################
################################################################################

#### Load Point Data

# Source: https://data.cityofnewyork.us/Education/School-Point-Locations/jfju-ynrr

(loc <- path.expand(paste0(getwd(), "/meta_data/Public_School_Locations/")))
file.exists( paste0(loc, 'Public_Schools_Points_2011-2012A.shp') )

shp <- readOGR(loc, layer="Public_Schools_Points_2011-2012A")
(pstr <- proj4string(shp))

summary(shp)

pj <- data.frame(project(shp@coords, pstr, inverse=TRUE))
pj <- cbind(pj, shp@data)

pjr <- select(pj, lat=X2, lng=X1, dbn=ATS_CODE, points_schoolname=SCHOOLNAME)
pjr$dbn <- stringr::str_trim(as.character(pjr$dbn))

# View(shp@data)
# View(shp@coords)

# plot(shp)

################################################################################
################################################################################
################################################################################

#### Load the Census data

# Source: This was downloaded with a custom ACS script not included with this package

# Join the Census Redux file
cen <- readRDS(paste0(prjdir, "census/", report_statabbr, "_ACS5_df_", 2016, ".Rds"))

# Do census tract numbers repeat across counties? Yes
ny=cen$tract[cen$county=="061"]       # Manhattan
queens=cen$tract[cen$county=="081"]   # Queens
bronx=cen$tract[cen$county=="005"]    # Bronx
kings=cen$tract[cen$county=="047"]    # Kings/Brooklyn
richmond=cen$tract[cen$county=="085"] # Richmond/Staten Island
nass=cen$tract[cen$county=="059"]     # Nassau

# Should be unique within county
all(sum(sum(duplicated(ny)), sum(duplicated(queens)), sum(duplicated(kings)), sum(duplicated(bronx)), sum(duplicated(richmond)))==0)

# Duplicated tract numbers have different data FYI. They are not repeated values across counties
sum(duplicated(c(ny, queens, kings,bronx,richmond)))



################################################################################
################################################################################
################################################################################

# Load DOE school file to save census tract and related school-level info

# Source: DOE

doe <- readxl::read_excel("meta_data/LCGMS_SchoolData_20180107_1420.xlsx")
doe <- filter(doe, !is.na(`ATS System Code`))

table(doe$City)
doe$City1 <- gsub("CITY", "", toupper(doe$City))
table(doe$City1)

# Construct state and county ids
doe$stateid <- 36
doe$countyid <- NA
doe$countyid <- with(doe,
                  ifelse(grepl("NEW YORK|NEW YORK|MANHATTAN", City1), "061",
                         ifelse(grepl("SYOSSET", City1), "059",
                                ifelse(grepl("BROOKLYN", City1), "047",
                                       ifelse(grepl("BRONX", City1), "005",
                                       ifelse(grepl("STATEN", City1), "085",
                                              ifelse(grepl("QUEENS|JACKSON|BAYSIDE|JAMAICA|ROCKAWAY|LONG ISLAND", City1), "081", NA)))))))


# Fix the tract for Syosset
doe$`Census Tract`[doe$City1=="SYOSSET"] <- 518800

# Every school has a county code
table(doe$countyid, useNA='always')

doe$`Census Tract` <- sprintf('%06d', as.numeric(doe$`Census Tract`))


# Build the GEOID field
doe$GEOID <- with(doe, paste0(stateid, countyid, `Census Tract`))

# Save a subset of fields
doerdx <- select(doe,
                 GEOID,
                 dbn=`ATS System Code`,
                 doe_name=`Location Name`,
                 doe_managed_by=`Managed By Name`,
                 doe_school_type=`Location Type Description`,
                 doe_grade_cat=`Location Category Description`,
                 doe_grade_levels=`Grades`,
                 #doe_grade_levels_final=`Grades Final`,
                 doe_open_year=`Open Date`,
                 doe_address=`Primary Address`,
                 doe_zip=Zip,
                 doe_city=City1,
                 doe_tract=`Census Tract`,
                 doe_bor_block=`Borough Block Lot`,
                 doe_superintendent_code=`Superintendent Location Code`)

table(doerdx$doe_city)

# Join Points data to DOE school roster

outersect <- function(x, y) {
    sort(c(setdiff(x, y),
           setdiff(y, x)))
}
# outersect(doerdx$dbn, pjr$dbn)


################################################################################
################################################################################
################################################################################

# Join the data


# DBNs in the school metafile
table(surveys_2016$dbn %in% doerdx$dbn)

# Missing appear to be mostly magnet schools (67):
# View(surveys_2016[!surveys_2016$dbn %in% doerdx$dbn, ])

# 2017
table(surveys_2017$dbn %in% doerdx$dbn)


## Points data matches for the survey file
table(surveys_2016$dbn %in% pjr$dbn) # Missing 175
table(surveys_2017$dbn %in% pjr$dbn) # Missing 176

# Missing Census tracts from DOE tract list:
table(doerdx$GEOID %in% cen$GEOID) # None

# Final merges:

# Join Census data to DOE list
sch_dat <- left_join(doerdx, cen)
nrow(sch_dat)==nrow(doerdx)

# Join the points data where available
sch_dat <- left_join(sch_dat, pjr)

# Join metadata to survey data
s16 <- right_join(sch_dat, surveys_2016)
nrow(s16)==nrow(surveys_2016)
saveRDS(s16, paste0(2016, '/data/merged_', 2016, '_all_vars.Rds'))
save(list="s16", file=paste0(2016, '/data/merged_', 2016, '_all_vars.RData'))

devtools::use_data(s16, overwrite = TRUE)


s17 <- right_join(sch_dat, surveys_2017)
nrow(s17)==nrow(surveys_2017)
saveRDS(s17, paste0(2017, '/data/merged_', 2017, '_all_vars.Rds'))
save(list="s17", file=paste0(2017, '/data/merged_', 2017, '_all_vars.RData'))

# Save out the master files to the package
devtools::use_data(s17, overwrite = TRUE)


# Function to calculate scores for a stub
df=s17
stub="s_q10a"

score_calc <- function(df, stub, percent_positive=TRUE) {
    names <- names(df)[grepl(stub, names(df))]
    cat( "Return formula: 100 * (", names[3], " + ", names[4], ") / (", names[1], " + ", names[2], " + ", names[3], " + ", names[4], ")\n", sep="")
    num <- apply( df[, c(names[3], names[4])], 1, sum, na.rm=TRUE)
    denom <- apply( df[,c(names[1:4])], 1, sum, na.rm=TRUE)
    denom[denom==0] <- NA
    return(100 * ( num / denom ))
}


stub="s_q5c"
df$test <- score_calc(df,  stub)
ttt=select(df, starts_with(stub), test)

# Missing Census Data by city
group_by(s17, doe_city) %>% summarize(n=n(),missing=sum(is.na(pct_pov)))

table(is.na(s17$pct_pov))
