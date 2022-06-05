* This is a new release, but was earlier submitted.  
* The license file has now been updated to comply with CRAN.  
* Some general references are now provided in the description.  
The website contains the most amount of information on this new methodology

There is no current publications, but I have 

This package covers an important development in new methodology I would like to
easily share with collaborators.



The package is tested using Github actions. The only remaining notes are about
GNU being needed with Rcpp and the size of the libraries, which seems to be
needed to make it work. Notes are summarized below.

## Github actions results

macOS-latest (release)

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors | 0 warnings | 1 note 


windows-latest (release)

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors | 0 warnings | 1 note

ubuntu-latest (release)

❯ checking installed package size ... NOTE
    installed size is 23.2Mb
    sub-directories of 1Mb or more:
      libs  22.7Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors | 0 warnings | 2 notes 