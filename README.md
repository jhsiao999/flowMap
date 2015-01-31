# flowMap

'flowMap' offers a distance metric for mapping cell populations in flow cytometry data. 
Our metric is distribution-free and can accommodate large number of feature markers. 
In a typical flow cytometry workflow, cell population mapping is a key process after
identifying cell populations in each flow cytometry samples to establish equivalence
of the cell populations across samples or to quantify biological similarity between 
cell types. 'flowMap' is a stand-alone method that can be easily incorporated with
any flow cytometry workflow.


## Preparing data for 'flowMap'

* Input: a set of gated sample that are stored in txt format. Each sample data file
is in a table format consisting of marker expression levels (rows) by marker channels
(columns). The last column of the file is required to be the cell population identifying
numbers generated from the gating method. 

* Output: a dissimilarity matrix comparing cell populations across samples. 


## Updating to the latest version of flowMap

You can track (and contribute to) the development of 'flowMap' at https://github.com/JoyceHsiao/flowMap. 

1. To install the release version of 'flowMap':
   ```R
   source("http://bioconductor.org/biocLite.R")
   biocLite("flowMap")
   ```

2. To install the development version of 'flowMap':
   ```R
   library("devtools")
   install_github("joycehsiao/flowMap")
   
   # If you don't have 'devtools' installed already
   install.packages("devtools")
   library("devtools")
   install_github("joycehsiao/flowMap")
   ```

