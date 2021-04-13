#' Retrieves the capturing dates of an image
#' 
#' Gets the date field from a file name
#' 
#' @param fnames a \code{character} vector of names
#' @param format a \code{character} specifying the date format to be returned
#' 
#' @returns a vector of \code{Date}s

get_dates <- function(fnames, format = "%Y%j"){
  as.Date(gsub(".*?([0-9]{1,7}).*$", "\\1", fnames), format)
}


#' Retrieves the capturing dates of an image
#' 
#' Gets the date field from a layer name
#' 
#' @param r a \code{Raster*} with layer names
#' @param format a \code{character} specifying the date format to be returned
#' 
#' @returns a vector of \code{Date}s

get_dates_from_layer <- function(r, format = "%Y%j"){
  as.Date(gsub(".*?([0-9]{1,7}).*?", "\\1", names(r)), format)
}
