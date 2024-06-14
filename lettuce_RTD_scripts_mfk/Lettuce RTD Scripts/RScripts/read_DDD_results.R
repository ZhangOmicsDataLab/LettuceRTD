#' Read a comma delimited files in 3D RNA-seq result folder
#'
#' @description 
#' This function reads comma delimited files (*.csv) in 3D RNA-seq result folder directly from zipped 3D RNA-seq output file.
#'  
#' @param zip.path A path to a 3D RNA-seq output file
#' @param job.id A 3D RNA-seq job id that outputs are saved in.
#' @param file.name A csv file name in result folder

read.DDD.results <- function(zip.path, job.id, file.name) {
  if (file.exists(paste0({{zip.path}})) == TRUE) {
    read.table(
      unz(description = paste0({{ zip.path }}),
          filename = paste0("srv/3drnaseq/", {{ job.id }}, "/result/", {{ file.name }})),
      header = T,
      quote = "\"",
      sep = ","
    )
  }
  else {
    print(paste("ERROR:", basename({{zip.path}}), "file does NOT exist in the directory!"))
    print(paste("Please make sure you typed the zip file's path or name correct."))
  }
}