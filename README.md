# CThripsos an R package to analyse chromothripsis in single cells. 

CThripsos is an easy to use R package that allows to analyse scDNA-seq data and to provide different scores that indicate the presence of chromotripsis in single cell CNV profiles.

# Installation
`devtools::install_github("gonzaparra/CThripsos")`

# Dependencies
ggplot2

# Minimum code
`library(CThripsos)`

`Segments_Matrix <- read.table("file.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)`

`CThripsosObject<-CreateCThripsosObject(Segments_Matrix)`

`window_length <-50000000`

`min_cnv_changes=10`

`min_consec_cnvs <- 1`

`CThripsosObject<-Calculate_CT_Cells(CThripsosObject, window_length, min_cnv_changes, min_consec_cnvs)`

# You can find an example of how to use the package at:
https://github.com/gonzaparra/CThripsos/tree/master/Examples
