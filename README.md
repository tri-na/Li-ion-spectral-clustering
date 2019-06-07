# Li-ion-spectral-clustering
Spectral clustering of mxrd spectrum for Li

Software: R version 3.3.3 
Package required: kernlab

### Installation guide
Follow the standard procedure at CRAN website of installing R and related packages

### Instructions for use and expected output
For spectral clustering, the function of "specc" from R library 'kernlab' is utilzied to cluster the different Li-ion compounds based on the mxrd spectrum. This code divisively as opposed to agglomerate in hierarhical clustering and iteratively split the larger cluster of Li-ion compounds into 2 until we obtained 8 clusters in total.

#### Input: Liion_comp_528.csv
#### expected output
    spectral clustering 
### expected runtime
<10 min on 8GM RAM, Intel i5-6300U CPU@2.4GHz
