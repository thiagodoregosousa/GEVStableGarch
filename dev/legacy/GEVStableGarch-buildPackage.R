
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  Several                 Commands to create the package to send to CRAN
################################################################################




--------------------------------------------------------------
# STEP: Clean the Environment
rm(list=ls()) 


--------------------------------------------------------------
# STEP: Declare directory that contain the files
my.files.directory = '/Users/thiago/dropbox/Pappers/GEVStableGarch/GEVStableGarch/'


--------------------------------------------------------------
# STEP: Declare directory that will contain the package
cran.directory = '/Users/thiago/dropbox/Pappers/GEVStableGarch/GEVStableGarch/CRAN_versions/'


--------------------------------------------------------------
# STEP: Version of package
my.version = '1.1'



--------------------------------------------------------------
# STEP: Create a folder with this name at
# '/Users/thiago/dropbox/Pappers/GEVStableGarch/CRAN_versions/'

  
--------------------------------------------------------------
# STEP: Create the Package Structure (man pages, function files, Description, etc)
# Declare list of files that will be contained inside the package
# If you include/delete a new file in the package do not forget to update this list
my.list.of.files = c("class-GEVSTABLEGARCH.R",
                     "class-GEVSTABLEGARCHSPEC.R",
                     "dist-skstd.R",
                     "dist-gat.R",
                     "GEVStableGarch-armaGarchDist.R",
                     "GEVStableGarch-filter.R",
                     "GEVStableGarch-fit.R",
                     "GEVStableGarch-getFormula.R",
                     "GEVStableGarch-getOrder.R",
                     "GEVStableGarch-getStart.R",
                     "GEVStableGarch-onAttach.R",
                     "GEVStableGarch-package.R",
                     "GEVStableGarch-select.R",
                     "GEVStableGarch-sim.R",
                     "GEVStableGarch-spec.R",
                     "GEVStableGarch-stationarityAparch.R",
                     "methods-show.R")

my.list.of.files.complete.adress = paste(my.files.directory, my.list.of.files, sep = '')

# DO NOT RUN IT TWICE BECAUSE YOU CAN LOOSE YOUR WORK 
#    package.skeleton(name='GEVStableGarch',
#                 path = paste(cran.directory,my.version,sep = ''),
#                 code_files = my.list.of.files.complete.adress,
#                 force = FALSE)



--------------------------------------------------------------
# STEP: Make necessary modification on files
  # Put DEBUG = FALSE in gsFit function and take it off
  # from the function parameters, making it equal to FALSE in the beggining of the function.
  # change the name of functions inside file armaGarchDist from "stable::dstable.quick"
  # to "GSgarch.dstable" (see file GEVStableGarch-onAttach )

  
  
  
--------------------------------------------------------------
# STEP: Fill in the package documentation
  # Help pages ( .Rd files)
  # DESCRIPTION
  # NAMESPACE
  # Put new references in the inst/doc/bibliography.txt file first and then copy to the
  # desired file

  
  
--------------------------------------------------------------
# STEP: Create the tar.gz file by executing the two commands on the terminal
#  < cd /Users/thiago/dropbox/pappers/GEVStableGarch/GEVStableGarch/CRAN_versions/1.1 > 
#  < R --vanilla CMD build GEVStableGarch > 
  
  
--------------------------------------------------------------  
# STEP: Check your package by running the command
#  Install the latest 'stable' R version since the CRAN
#  maintainers will use it to check your package
#  1: < R --vanilla CMD check GEVStableGarch_1.1.tar.gz > 
#  2: < R CMD check --as-cran GEVStableGarch_1.1.tar.gz  > 
#  3: Test using the current R-devel version and R-release
#  versions using the website:
#  http://win-builder.r-project.org/upload.aspx

  
  
--------------------------------------------------------------  
# STEP: Copy additional files to the package
# changeLog to the root
  
 
  
--------------------------------------------------------------  
# STEP: Final test before sending to CRAN
# Install the package on both mac and windows and run 
# the GEVStableGarch-TestAfterBuildPackage.R file inside R (not R studio)



--------------------------------------------------------------  
# STEP: Send to CRAN: 
# https://cran.r-project.org/web/packages/policies.html#Submission
# Notes: 
# 1 - If there is an warning you cannot eliminate, explain the
# reason for that on the submission form.
# 2 - Before submitting a package update, consult the CRAN 
# check the error report page at 
# https://cran.r-project.org/web/checks/check_results_GEVStableGarch.html
# 3 - When emailing cran maintainers always send a copy (cc) to 
# "CRAN@R-project.org" <CRAN@r-project.org>
# 4 - Notes we are unable to solve
#   4.1 - Package 'stable' is only available at http://www.robustanalysis.com and
#   therefore cannot be checked. Most users of package GEVStableGarch own the 'stable' 
#   library and hence, the enhances sections is really important.

  
--------------------------------------------------------------  
# ADIVICES FOR CORRECTING SOME ERRORS
--------------------------------------------------------------

# CHECKING FOR NON-ASCII CHARACTERS: Use function showNonASCIIfile(file) from package tools
library(tools)
file.to.find.non.ascii = "/GEVStableGarch/man/gsFit.Rd"
showNonASCIIfile(paste(cran.directory, my.version, file.to.find.non.ascii,sep = ''))

# CHECK THE PACKAGE EVERYTIME YOU PUT NEW LATEX CODE INTO IT. WHY ? 
# The R CMD check does not give you advice about where is the error.

# DO NOT PUT BLANK LINES BETWEEN \eqn OR \deqn COMMANDS

# BUID PACKAGE TO CORRECT ERRORS

# DO NOT USE THE \text or \texttt MACRO INSIDE THE \eqn ENVIRONMENT




