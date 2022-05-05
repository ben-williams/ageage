#' Setup folder structure
#'
#' Creates a common folder structure for ADMB assessments
#'
#' @param year
#' #'
#' @return
#' @export setup
#'
#' @examples
#' \dontrun{
#' setup(2020)
#'}
setup <- function(year){
  
  if(!dir.exists(here::here(year, "data", "output"))) {
    dir.create(here::here(year, "data", "output"), recursive = TRUE)
  }
  
    dirs = "models"
    folders = "ageage"
    
    for(i in 1:length(dirs)){
      if(dir.exists(here::here(year, "data", dirs[i])) == FALSE){
        dir.create(here::here(year, "data", dirs[i]), recursive = TRUE)
      }
    }
    for(i in 1:length(folders)){
      dir.create(here::here(year, "data", "models", folders[i]))
    }
    
    file.copy(system.file("tpl", "ageage.tpl", package = "ageage"),
              here::here(year, "data", "models", "ageage"))
    
  #   file.copy(system.file("models", "allometric.tpl", package = "groundfishr"),
  #             here::here(year, "data", "models", "allometric"))
  #   
  #   file.copy(system.file("models", "vbl.tpl", package = "groundfishr"),
  #             here::here(year, "data", "models", "vonb"))
  #   
  #   file.copy(system.file("models", "wvbl.tpl", package = "groundfishr"),
  #             here::here(year, "data", "models", "wvonb"))
  #   
  #   file.copy(system.file("models", "wvbl.ctl", package = "groundfishr"),
  #             here::here(year, "data", "models", "wvonb"))
  #   
  #   file.copy(system.file("models", "lengthsd.tpl", package = "groundfishr"),
  #             here::here(year, "data", "models", "length_sd"))
  # }
  # dir.create(here::here(year, "proj"))
  # dir.create(here::here(year, "apportionment"))
  
}