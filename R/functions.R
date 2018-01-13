#' Calculate percent positive response for a question in the New York City School Survey data
#' @param df data.frame; s16 or s17
#' @param stub character; The question stub of interest, e.g. "s_q5c", which is the stub for "5c. My teachers will always listen to students' ideas."
#' @return A numeric vector of proportions with console output indicating the selected variables
#' @examples
#' \dontrun{
#' df <- s17
#' df$test <- score_calc(df,  stub="s_q5c")
#' }
#' @export

score_calc <- function(df, stub) {
    names <- names(df)[grepl(stub, names(df))]
    message( "Return formula: 100 * (", names[3], " + ", names[4], ") / (", names[1], " + ", names[2], " + ", names[3], " + ", names[4], ")\n", sep="")
    num <- apply( df[, c(names[3], names[4])], 1, sum, na.rm=TRUE)
    denom <- apply( df[,c(names[1:4])], 1, sum, na.rm=TRUE)
    denom[denom==0] <- NA
    return(100 * ( num / denom ))
}
