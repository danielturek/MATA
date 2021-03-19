library(devtools)
library(roxygen2)

baseDir <- '~/github/MATA/'
setwd(baseDir)

document(paste0(baseDir, 'MATA'))
document(paste0(baseDir, 'MATA'))
system(paste0('R CMD BUILD ', baseDir, 'MATA'))

check(paste0(baseDir, 'MATA'))

suppressMessages(try(remove.packages('MATA'), silent = TRUE))
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
message('installing package version ', gsub('\\.tar\\.gz$', '', lastTarFile))
system(paste0('R CMD install ', lastTarFile))

q('no')

## now restart R

library(MATA)

?MATA
?mata
?mata.wald
?dmata.wald
?pmata.wald


mata.wald
dmata.wald
pmata.wald








