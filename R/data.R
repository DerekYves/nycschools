#' 2017 NYC School Survey Results Joined with Census and School Meta Data
#'
#' This is a subset of parent and student survey questions from the New York City
#' School Survey. It joins a subset of Census data and information about each
#' participating school from a meta file provided by the DOE. See the git repository
#' for this project for additional source files and the data dictionary:
#' https://github.com/DerekYves/nycschools
#'
#' @source  http://schools.nyc.gov/Accountability/tools/survey/2017Survey.htm
#'
#'
#' @format A data frame with 1,137 rows and 475 variables:
#' \describe{
#'   \item{p_q*}{Responses from the Parent survey}
#'   \item{s_q*}{Responses from the Student survey}
#'   \item{colab}{Composite score for the school: Collaborative Teachers}
#'   \item{eff_leaders}{Composite score for the school: Effective leadership}
#'   \item{rig_instr}{Composite score for the school: Rigorous Instruction}
#'   \item{supp_env}{Composite score for the school: Supportive Environment}
#'   \item{stg_ties}{Composite score for the school: Strong Family/Comminity Ties}
#'   \item{trust}{Composite score for the school: Trust Score}
#'   \item{doe*}{Variables from the DOE school meta file}
#'   \item{pct_unv_deg}{Census: Percent of school Census tract with a university degree}
#'   \item{pct_under_hs}{Census: Percent of school Census tract with less than HS education}
#'   \item{pct_pov}{Census: Percent of school Census tract below poverty line}
#'   \item{pct_inc_100k_up}{Census: Percent of school Census tract with income 100k or higher}
#'   \item{pct_inc_150k_up}{Census: Percent of school Census tract with income 150k or higher}
#'   \item{est_med_hinc}{Census: Estimated median household income of school Census tract}
#'   \item{ratio_med_hinc}{Census: Ratio of median household income within tract to all school tracts in sample}
#'   \item{pct_any_hc_u65}{Census: Percent of school Census tract with with public or private healthcare}
#'   \item{pct_prv_hc_u65}{Census: Percent of school Census tract with with private healthcare}
#'   \item{pct_black/white/amind/asian/nhpi/oth/latin/pct_mlt_r/nonwht}{Census: Percent of school Census tract within by Census race/rthnicity categories}
#'   \item{sum_*}{Census: Raw counts for a given category within the school's Census tract}
#'   \item{lat}{DOE Points File: Latitude of school}
#'   \item{lng}{DOE Points File: Longitude of school}
#'   \item{year}{Year of Survey}
#'   ...
#' }
"s17"

#' 2016 NYC School Survey Results Joined with Census and School Meta Data
#'
#' This is a subset of parent and student survey questions from the New York City
#' School Survey. It joins a subset of Census data and information about each
#' participating school from a meta file provided by the DOE. See the git repository
#' for this project for additional source files and the data dictionary:
#' https://github.com/DerekYves/nycschools
#'
#' @source  http://schools.nyc.gov/Accountability/tools/survey/2016Survey.htm
#'
#'
#' @format A data frame with 1,118 rows and 463 variables:
#' \describe{
#'   \item{p_q*}{Responses from the Parent survey}
#'   \item{s_q*}{Responses from the Student survey}
#'   \item{par_sat}{Composite score for the school: Parent Satisfaction}
#'   \item{stu_sat}{Composite score for the school: Student Satisfaction}
#'   \item{tch_sat}{Composite score for the school: Teacher Satisfaction}
#'   \item{colab}{Composite score for the school: Collaborative Teachers}
#'   \item{eff_leaders}{Composite score for the school: Effective leadership}
#'   \item{rig_instr}{Composite score for the school: Rigorous Instruction}
#'   \item{supp_env}{Composite score for the school: Supportive Environment}
#'   \item{stg_ties}{Composite score for the school: Strong Family/Comminity Ties}
#'   \item{trust}{Composite score for the school: Trust Score}
#'   \item{doe*}{Variables from the DOE school meta file}
#'   \item{pct_unv_deg}{Census: Percent of school Census tract with a university degree}
#'   \item{pct_under_hs}{Census: Percent of school Census tract with less than HS education}
#'   \item{pct_pov}{Census: Percent of school Census tract below poverty line}
#'   \item{pct_inc_100k_up}{Census: Percent of school Census tract with income 100k or higher}
#'   \item{pct_inc_150k_up}{Census: Percent of school Census tract with income 150k or higher}
#'   \item{est_med_hinc}{Census: Estimated median household income of school Census tract}
#'   \item{ratio_med_hinc}{Census: Ratio of median household income within tract to all school tracts in sample}
#'   \item{pct_any_hc_u65}{Census: Percent of school Census tract with with public or private healthcare}
#'   \item{pct_prv_hc_u65}{Census: Percent of school Census tract with with private healthcare}
#'   \item{pct_black/white/amind/asian/nhpi/oth/latin/pct_mlt_r/nonwht}{Census: Percent of school Census tract within by Census race/rthnicity categories}
#'   \item{sum_*}{Census: Raw counts for a given category within the school's Census tract}
#'   \item{lat}{DOE Points File: Latitude of school}
#'   \item{lng}{DOE Points File: Longitude of school}
#'   \item{year}{Year of Survey}
#'   ...
#' }
"s16"

