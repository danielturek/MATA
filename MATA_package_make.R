
library(devtools)
##setwd('/Users/dturek/Dropbox/Berkeley/R/2_MATA_package')
setwd('~/github/MATA')

descList <- list(
    Title = 'Model-Averaged Tail Area Wald (MATA-Wald) Confidence Interval',
    Version = '0.4',
    `Authors@R` = "person('Daniel', 'Turek', email='danielturek@gmail.com', role=c('aut','cre'))",
    Description = 'Calculates Model-Averaged Tail Area Wald (MATA-Wald) confidence intervals, which are constructed using single-model estimators and model weights.',
    License = 'GPL-2'
)

## removes existing package
if('MATA' %in% dir()) unlink('MATA', recursive=TRUE, force=TRUE)

remove.packages('MATA')
usethis::create_package('MATA', fields=descList, rstudio=FALSE)

setwd('~/github/MATA')
system('cp MATA_package_source.R MATA/R/')

document('MATA')
build('MATA')
check('MATA')
q('no')

## shell: $ R CMD INSTALL MATA_0.1.tar.gz
## shell: $ R CMD INSTALL MATA_0.2.tar.gz
## shell: $ R CMD INSTALL MATA_0.3.tar.gz
## shell: $ R CMD INSTALL MATA_0.4.tar.gz

##setwd('/Users/dturek/Dropbox/Berkeley/R/2_MATA_package')
setwd('~/github/MATA')
library('MATA')

?mata.wald







