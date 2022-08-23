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

# If you find this package useful please cite:
"Single cell multi-omics analysis of chromothriptic medulloblastoma highlights genomic and transcriptomic consequences of genome instability". R Gonzalo Parra, Moritz J Przybilla, Milena Simovic, Hana Susak, Manasi Ratnaparkhe, John KL Wong, Verena Körber, Philipp Mallm, Martin Sill, Thorsten Kolb, Rithu Kumar, Nicola Casiraghi, David R Norali Ghasemi, Kendra Korinna Maaß, Kristian W Pajtler, Anna Jauch, Andrey Korshunov, Thomas Höfer, Marc Zapatka, Stefan M Pfister, Oliver Stegle, Aurélie Ernst.

biorXiv 2021. doi: https://doi.org/10.1101/2021.06.25.449944
