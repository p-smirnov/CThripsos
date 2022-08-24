library(CThripsos)

Segments_Matrix <- read.table("/Users/parra/Desktop/StegleGroup/Meduloblastoma/cnv.cells.mat.20kb-bin.corrected.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
dim(Segments_Matrix)

# Cluster annotation
clusters <-read.table("/Users/parra/Desktop/StegleGroup/Meduloblastoma/STP-Nuclei.txt", sep = "\t", stringsAsFactors = F, header = T)
clusters[,"cell_id"]<-paste0("X",clusters[,"cell_id"])

# We subselect only the cells that are present in the annotated clusters
Segments_Matrix<-Segments_Matrix[,clusters[,"cell_id"]]
dim(Segments_Matrix)

Segments_Matrix<- SortMatrixChrsNumerically(Segments_Matrix)

# We create metadata objects and store the CNV matrix in a CThripsos object
CThripsosObject<-CreateCThripsosObject(Segments_Matrix)
# rm(Segments_Matrix)

window_length <-50000000
min_cnv_changes=10
min_consec_cnvs <- 1

CThripsosObject<-Calculate_CT_Cells(CThripsosObject, window_length, min_cnv_changes, min_consec_cnvs)

CellsClusters<-cbind(colnames(Segments_Matrix), 1:length(colnames(Segments_Matrix)))
colnames(CellsClusters)<-c("cell_id", "cluster")
# We calculate CT for all metacells that were stored in the CThripsosObject

# We create metacells based on the clusters that were provided in the earlier step
CThripsosObject<-CreateMetacells(CThripsosObject, CellsClusters)

CThripsosObject<-Calculate_CT_Metacells(CThripsosObject, window_length, min_cnv_changes, min_consec_cnvs)

plot_MetacellsCT(CThripsosObject)




