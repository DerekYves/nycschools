
# Score calculation
score_calc <- function(df, stub, percent_positive=TRUE) {
    names <- names(df)[grepl(stub, names(df))]
    cat( "Return formula: 100 * (", names[3], " + ", names[4], ") / (", names[1], " + ", names[2], " + ", names[3], " + ", names[4], ")\n", sep="")
    num <- apply( df[, c(names[3], names[4])], 1, sum, na.rm=TRUE)
    denom <- apply( df[,c(names[1:4])], 1, sum, na.rm=TRUE)
    denom[denom==0] <- NA
    return(100 * ( num / denom ))
}
