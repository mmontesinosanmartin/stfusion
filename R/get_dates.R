#' Retrieves the capturing dates of an image
#' 
#' Gets the date field from a file name
#' 
#' @param filenames a \code{character} vector of file names
#' @param format a \code{character} specifying the date format
#' 
#' @returns a vector of \code{Date}s

get_dates <- function(filenames, format){
  as.Date(gsub(".*?([0-9]{1,7}).*$", "\\1", filenames), format)
}
