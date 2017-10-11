RESUBMISSION, updated to version 0.1.2 from 0.1.1

Processed the following remarks kindly provided by Swetlana Herbrandt <herbrandt@statistik.tu-dortmund.de>

## Remark doi
Thanks, please omit the space after 'doi:':
<doi:10.1093/bioinformatics/btu091>
* removed spaces

## Remark writing to home dir
Your examples
map_export(MF.obj = MF.obj, file = "consensus")
still writes in the user's home directory.

* Added tempdir to dontrun example

## Other updates
* Changed 'Create Consensus Genetic..' in Description field of DESCRIPTION to 'Construct Consensus Genetic..'
* Update map_export parameter documentation

## R CMD check --as-cran results on mapfuser_0.1.2.tar.gz
There were no ERRORs or WARNINGs or NOTEs 
Maintainer:  Dennis van Muijen <d.van.muijen@rijkzwaan.nl> 

## Test environments
* local Ubuntu 14.04 install, R 3.4.1
* local Ubuntu 14.04 install, R-devel
* win-builder (R-devel)
* local windows 7, R 3.4.1


##### Resubmission 2
RESUBMISSION, updated to version 0.1.2 from 0.1.1

Processed the following remarks kindly provided by Swetlana Herbrandt <herbrandt@statistik.tu-dortmund.de>

## Remark doi
Thanks, please omit the space after 'doi:':
<doi:10.1093/bioinformatics/btu091>

* removed spaces

## Remarkt file write to home directory
Your examples
map_export(MF.obj = MF.obj, file = "consensus")
still writes in the user's home directory.

* Added tempdir to dontrun example

## Other updates
* Changed 'Create Consensus Genetic..' in Description field of DESCRIPTION to 'Construct Consensus Genetic..' 
* Changed 'genomic' to 'physical distance' in Description field of DESCRIPTION
* Update map_export parameter documentation
* Removed commented out code from vignette

## R CMD check --as-cran results on mapfuser_0.1.2.tar.gz
There were no ERRORs or WARNINGs or NOTEs other than a NOTE related to package maintainer:
NOTE
Maintainer:  Dennis van Muijen <d.van.muijen@rijkzwaan.nl> 

## Test environments
* local Ubuntu 14.04 install, R 3.4.1
* local Ubuntu 14.04 install, R-devel
* win-builder (R-devel)
* local windows 7, R 3.4.1