#' Download the BioGRID COVID-19 dataset
#' @example{\dontrun{
#'
#'try({
#'  usethis::use_zip(
#'    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-PROJECT-covid19_coronavirus_project-LATEST#'.zip",
#'    destdir="data-raw"
#'  )
#'})
#'
#'unzip(
#'  zipfile="data-raw/BIOGRID-PROJECT-covid19_coronavirus_project-LATEST.zip",
#'  exdir="data-raw/BIOGRID-PROJECT-covid19_coronavirus_project-LATEST"
#')
#'
#' }}
