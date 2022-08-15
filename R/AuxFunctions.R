
getAnnotations<-function(gene_name)
{
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  listAttributes(ensembl, page="feature_page")
  
  annot <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                 filters = 'hgnc_symbol', 
                 values = gene_name, 
                 mart = ensembl)
  
  annot <- annot[which(annot[,"chromosome_name"] %in% c(1:22,"X","Y")),]
  
  unique_gene_ids<-which(!duplicated(annot[,"hgnc_symbol"]))
  annot <-annot[unique_gene_ids, ]
  rownames(annot) <-annot[,"hgnc_symbol"]
  
  return(annot)
}

getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg19"){
  g <- getBSgenome(genome, masked=FALSE)
  data.frame(chrom=1:24, chrXY=c(1:22,'X','Y'), length=seqlengths(g)[1:24])
}

# get cyto bands
getCytoBand <- function( genome="hg19" ){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- ucscTableQuery(mySession, table="cytoBandIdeo") #cytoBand
  tbl <- getTable(obj)
  .chrAsNum(tbl)
}


GeneCNVProfile<-function (gene_chr, gene_start, gene_end, segment_coords, Segments_Matrix, title)
{
  # # example
  # gene_chr=2
  # gene_start=121493199
  # gene_end=121750229
  
  
  # take the segments within the gene is contained
  gene_segments<-c(which(as.numeric(segment_coords[,1]) <= gene_start & as.numeric(segment_coords[,2]) > gene_start & segment_chromosomes==gene_chr), which(as.numeric(segment_coords[,1]) < gene_end & as.numeric(segment_coords[,2]) >= gene_end & segment_chromosomes==gene_chr))
  gene_segments<-unique(gene_segments)
  
  segment_coords[gene_segments,]
  # fields::image.plot(as.matrix(Segments_matrix8[gene_segments[1]:gene_segments[-1],]), xlab="chromosomes", ylab="cells", axes=F)
  if(length(gene_segments)==1){gene_segments=c(gene_segments, gene_segments)}
  colMeans(t(Segments_Matrix[gene_segments[1]:gene_segments[-1],]))
  
  gene_cnv_profile<-as.vector(unlist(Segments_Matrix[gene_segments[1]:gene_segments[-1],]))
  gene_cnv_counts<-c()
  for(i in 0:range(gene_cnv_profile)[2])
  {
    gene_cnv_counts <-c(gene_cnv_counts, length(which(gene_cnv_profile==i)) / length(gene_cnv_profile))
  }
  
  colors<-rep("red", length(0:range(gene_cnv_profile)[2]))
  colors[which(0:range(gene_cnv_profile)[2]==2)]="gray"
  colors[which(0:range(gene_cnv_profile)[2]>2)]="blue"
  
  barplot(gene_cnv_counts, names.arg= 0:range(gene_cnv_profile)[2], col=colors, xlab = "Copy-number", ylab = "% of cells", main=title)
}

GeneCNVProfile_FISH<-function(fish_data, xrange, title)
{
  colors<-rep("red", length(fish_data))
  colors[which(xrange==2)]="gray"
  colors[which(xrange>2)]="blue"
  barplot(fish_data , names.arg= xrange, col=colors, xlab = "Copy-number", ylab = "% of cells", main = title)
}

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
