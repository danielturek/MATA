setwd('/Users/dturek/Dropbox/Berkeley/R/2_MATA_package')
library(devtools)
descList <- list(
    Title = 'Model-Averaged Tail Area Wald (MATA-Wald) Confidence Interval',
    Version = '0.3',
    `Authors@R` = "person('Daniel', 'Turek', email='dturek@maths.otago.ac.nz', role=c('aut','cre'))",
    Description = 'Calculates Model-Averaged Tail Area Wald (MATA-Wald) confidence intervals, which are constructed using single-model estimators and model weights.',
    License = 'GPL-2'
)
if('MATA' %in% dir()) unlink('MATA', recursive=TRUE, force=TRUE)  # removes existing package
remove.packages('MATA')
create('MATA', description=descList, rstudio=FALSE)
system('cp MATA_package_source.R MATA/R/')
document('MATA')
build('MATA')
check('MATA')
q('no')

## shell: $ R CMD INSTALL MATA_0.3.tar.gz

setwd('/Users/dturek/Dropbox/Berkeley/R/2_MATA_package')
library('MATA')
?mata.wald







