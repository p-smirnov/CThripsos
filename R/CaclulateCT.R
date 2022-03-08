# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

CT_scoring_single<-function(cell, window_length, min_cnv_changes, min_consec_cnvs, CThripsosObject)
{
  chromatriptic_cell<-FALSE
  chromosomes<- unique(CThripsosObject$Annotations$segment_chromosomes)

  chromatriptic_chromosomes<-c()
  bins_windows<-rep(0, nrow(CThripsosObject$CNVMatrix))

  for (chromosome_i in chromosomes)
  {
    chromatriptic_chromosome <-0
    chr_windows<-0
    segments_in_chromosome_i=which(CThripsosObject$Annotations$segment_chromosomes==chromosome_i)

    compressed_cell<-rle(cell[segments_in_chromosome_i])
    if (length(compressed_cell$values)<10) {
      # print("skipping chr...")
      chromatriptic_chromosomes=c(chromatriptic_chromosomes, chromatriptic_chromosome)
      next
    }

    for (segment_i in segments_in_chromosome_i)
    {
      window_limit <- segments_in_chromosome_i[which(as.numeric(CThripsosObject$Annotations$segment_coords[segments_in_chromosome_i,2]) - as.numeric(CThripsosObject$Annotations$segment_coords[segment_i, 1]) > window_length)[1]]

      if(!is.na(window_limit)) # this check is for the case we reached the border condition
      {
        window_limit_length=as.numeric(CThripsosObject$Annotations$segment_coords[window_limit,2]) - as.numeric(CThripsosObject$Annotations$segment_coords[segment_i, 1])
        # This is the segment where the changes need to be quantified.
        if(segment_i>window_limit)
        {
          print(paste("shit ",segment_i, " ",window_limit))
        }
        CNV_changes <-cell[segment_i:window_limit]

        # Take NAs out
        CNV_changes<-CNV_changes[which(!is.na(CNV_changes))]
        names(CNV_changes) <- NULL

        # We slide the window and count again until the end.
        compressed_cnv <- rle(CNV_changes)
        window_changes<- length(compressed_cnv$values)
        cnv_states<-unique(compressed_cnv$values)

        windows_to_delete<-which(compressed_cnv$lengths>=min_consec_cnvs)
        compressed_cnv$values<-compressed_cnv$values[windows_to_delete]
        compressed_cnv$lengths<-compressed_cnv$lengths[windows_to_delete]

        compressed_cnv<-inverse.rle(compressed_cnv)
        compressed_cnv<-rle(compressed_cnv)
        window_changes<- length(compressed_cnv$values)
        cnv_states<-unique(compressed_cnv$values)

        chr_windows <- chr_windows +1
        if(window_changes>= min_cnv_changes)
        {
          chromatriptic_cell=TRUE
          chromatriptic_chromosome=chromatriptic_chromosome+1
          bins_windows[segment_i:window_limit]<-bins_windows[segment_i:window_limit] +1
        }
      }
    }
    if(chromatriptic_chromosome >0)
    {
      chromatriptic_chromosome <- chromatriptic_chromosome/chr_windows
      bins_windows[segments_in_chromosome_i]<-bins_windows[segments_in_chromosome_i]/chr_windows
    }
    chromatriptic_chromosomes=c(chromatriptic_chromosomes, chromatriptic_chromosome)
  }

  print(chromatriptic_chromosomes)
  return(list(chromatriptic_chromosomes=chromatriptic_chromosomes, bins_windows=bins_windows))
}


Calculate_CT_Metacells<-function(CThripsosObject, window_length, min_cnv_changes, min_consec_cnvs)
{
  CT_MetacellsChrs<-c()
  CT_MetacellsBins<-c()

  for(i in 1:nrow(CThripsosObject$Metacells$MetacellsMatrix))
  {
    print(paste("index", i, "..."))
    cell<-CThripsosObject$Metacells$MetacellsMatrix[i,]
    chromatriptic_chromosomes<-CT_scoring_single(cell, window_length, min_cnv_changes, min_consec_cnvs, CThripsosObject)
    CT_MetacellsChrs<-rbind(CT_MetacellsChrs, chromatriptic_chromosomes$chromatriptic_chromosomes)
    CT_MetacellsBins<-rbind(CT_MetacellsBins, chromatriptic_chromosomes$bins_windows)
  }

  print("finished calculating now storing objects..")
  rownames(CT_MetacellsChrs)<-rownames(CThripsosObject$Metacells$MetacellsMatrix)
  colnames(CT_MetacellsChrs)<-CThripsosObject$Annotations$chromosomes

  rownames(CT_MetacellsBins)<-rownames(CThripsosObject$Metacells$MetacellsMatrix)
  colnames(CT_MetacellsBins)<-colnames(CThripsosObject$Metacells$MetacellsMatrix)

  CThripsosObject$Metacells=c(CThripsosObject$Metacells, CT_MetacellsChrs=1)
  CThripsosObject$Metacells$CT_MetacellsChrs<-CT_MetacellsChrs

  CThripsosObject$Metacells=c(CThripsosObject$Metacells, CT_MetacellsBins=1)
  CThripsosObject$Metacells$CT_MetacellsBins<-CT_MetacellsBins

  CThripsosObject$Metacells=c(CThripsosObject$Metacells, window_length=1)
  CThripsosObject$Metacells$window_length<-window_length

  CThripsosObject$Metacells=c(CThripsosObject$Metacells, min_cnv_changes=1)
  CThripsosObject$Metacells$min_cnv_changes<-min_cnv_changes

  CThripsosObject$Metacells=c(CThripsosObject$Metacells, min_consec_cnvs=1)
  CThripsosObject$Metacells$min_consec_cnvs<-min_consec_cnvs

  return(CThripsosObject)
}


Calculate_CT_Cells<-function(CThripsosObject, window_length, min_cnv_changes, min_consec_cnvs)
{
  CT_CellsChrs<-c()
  CT_CellsBins<-c()

  for(i in 1:ncol(CThripsosObject$CNVMatrix))
  {
    print(paste("processing cell", i, "..."))
    cell<-CThripsosObject$CNVMatrix[,i]
    chromatriptic_chromosomes<-CT_scoring_single(cell, window_length, min_cnv_changes, min_consec_cnvs, CThripsosObject)
    CT_CellsChrs<-rbind(CT_CellsChrs, chromatriptic_chromosomes$chromatriptic_chromosomes)
    CT_CellsBins<-rbind(CT_CellsBins, chromatriptic_chromosomes$bins_windows)
  }

  rownames(CT_CellsChrs)<-colnames(CThripsosObject$CNVMatrix)
  colnames(CT_CellsChrs)<-CThripsosObject$Annotations$chromosomes

  rownames(CT_CellsBins)<-colnames(CThripsosObject$CNVMatrix)
  colnames(CT_CellsBins)<-rownames(CThripsosObject$CNVMatrix)

  CThripsosObject$CTData$CT_CellsChrs<-CT_CellsChrs
  CThripsosObject$CTData$CT_CellsBins<-CT_CellsBins

  CThripsosObject$CTData=c(CThripsosObject$CTData, window_length=1)
  CThripsosObject$CTData$window_length<-window_length

  CThripsosObject$CTData=c(CThripsosObject$CTData, min_cnv_changes=1)
  CThripsosObject$CTData$min_cnv_changes<-min_cnv_changes

  CThripsosObject$CTData=c(CThripsosObject$CTData, min_consec_cnvs=1)
  CThripsosObject$CTData$min_consec_cnvs<-min_consec_cnvs

  return(CThripsosObject)
}

# This function assumes that the columns contain bin annotation in the form chr:start-end
CreateCThripsosObject<-function(CNVMatrix)
{
  Annotations<-list()

  Annotations=c(Annotations, full_coords=1)
  Annotations$full_coords <- strsplit(rownames(CNVMatrix), split=":", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  Annotations$full_coords <- matrix(unlist(Annotations$full_coords), ncol=2, nrow=length(rownames(CNVMatrix)), byrow = T)

  Annotations=c(Annotations, segment_chromosomes=1)
  splitted_ids <- strsplit(rownames(CNVMatrix), split=":", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  splitted_ids <- matrix(unlist(splitted_ids), ncol=2, nrow=length(rownames(CNVMatrix)), byrow = T)
  Annotations$segment_chromosomes <- splitted_ids[,1]

  segment_coords <- Annotations$full_coords[,2]
  splitted_ids <- strsplit(segment_coords,  split="-", fixed = FALSE, perl = FALSE, useBytes = FALSE)

  Annotations=c(Annotations, segment_coords=1)
  Annotations$segment_coords <- matrix(unlist(splitted_ids), ncol=2, nrow=length(segment_coords), byrow = T)

  Annotations=c(Annotations, chromosomes=1)
  Annotations$chromosomes<- unique(Annotations$segment_chromosomes)

  CThripsosObject<-list()
  CThripsosObject=c(CThripsosObject, CNVMatrix=1)
  CThripsosObject$CNVMatrix<-CNVMatrix

  CThripsosObject=c(CThripsosObject, Annotations=1)
  CThripsosObject$Annotations<-Annotations

  CThripsosObject=c(CThripsosObject, CTData=1)
  CThripsosObject$CTData<-list()

  # This will be the matrix that will contain the CT data at the chr level.
  CT_CellsChrs<-matrix(NaN,nrow=ncol(CThripsosObject$CNVMatrix),ncol=length(CThripsosObject$Annotations$chromosomes))
  CThripsosObject$CTData=c(CThripsosObject$CTData, CT_CellsChrs=1)
  CThripsosObject$CTData$CT_CellsChrs<-CT_CellsChrs
  rm(CT_CellsChrs)
  rownames(CThripsosObject$CTData$CT_CellsChrs)<-colnames(CThripsosObject$CNVMatrix)
  colnames(CThripsosObject$CTData$CT_CellsChrs)<-CThripsosObject$Annotations$chromosomes

  CT_CellsBins<-matrix(NaN,nrow=nrow(CThripsosObject$CNVMatrix),ncol=ncol(CThripsosObject$CNVMatrix))
  CThripsosObject$CTData=c(CThripsosObject$CTData, CT_CellsBins=1)
  CThripsosObject$CTData$CT_CellsBins<-CT_CellsBins
  rm(CT_CellsBins)
  rownames(CThripsosObject$CTData$CT_CellsBins)<-rownames(CThripsosObject$CNVMatrix)
  colnames(CThripsosObject$CTData$CT_CellsBins)<-colnames(CThripsosObject$CNVMatrix)

  return(CThripsosObject)
}

# This function creates metacells based on cluster annotations
# col1 contains the cell ids, col2 contains cluster assignment
CreateMetacells<-function(CThripsosObject, clusters)
{
  Metacells<-c()

  for(cluster_i in unique(clusters[,"cluster"]))
  {
    print(paste0("creating metacell for cluster ",  cluster_i, "..."))
    metacell<- apply(as.matrix(CThripsosObject$CNVMatrix[,which(clusters[,"cluster"]==cluster_i)]), 1, function(x) as.numeric(names(which.max(table(as.numeric(x))))))
    Metacells<-rbind(Metacells, metacell)
  }
  rownames(Metacells)<-unique(clusters[,"cluster"])
  colnames(Metacells)<-rownames(CThripsosObject$CNVMatrix)

  CThripsosObject=c(CThripsosObject, Metacells=1)
  CThripsosObject$Metacells<-list()
  CThripsosObject$Metacells=c(CThripsosObject$Metacells, MetacellsMatrix=1)
  CThripsosObject$Metacells$MetacellsMatrix<-Metacells

  CThripsosObject$Metacells=c(CThripsosObject$Metacells, MetacellsClusters=1)
  CThripsosObject$Metacells$MetacellsClusters<-clusters

  print("Metacells object was added to CThripsosObject$Metacells...")

  return(CThripsosObject)
}

plot_MetacellsCT<-function (CThripsosObject, score_binary=T, rows=NULL)
{
  library(ggplot2)
  
  if(rows==NULL)
   {
    rows<-nrow(CThripsosObject$Metacells$MetacellsMatrix)
   }
  MetaCell_All_df <- c()
  plots <- list()
  p_i = 1
  for (ClonePlot in unique(CThripsosObject$Metacells$MetacellsClusters[,
                                                                       2])) {
    MetaCell_df <- as.data.frame(cbind(rep(ClonePlot, length(CThripsosObject$Annotations$segment_chromosomes)),
                                       colnames(CThripsosObject$Metacells$MetacellsMatrix)))
    colnames(MetaCell_df) <- c("Clone", "Chromosome")
    MetaCell_df$chrXY <- unlist(lapply(strsplit(MetaCell_df$Chromosome,
                                                ":"), "[[", 1))
    MetaCell_df$start <- unlist(lapply((strsplit(unlist(lapply(strsplit(MetaCell_df$Chromosome,
                                                                        "-"), "[[", 1)), ":")), "[[", 2))
    MetaCell_df$end <- unlist(lapply(strsplit(MetaCell_df$Chromosome,
                                              "-"), "[[", 2))
    MetaCell_df$CNV <- CThripsosObject$Metacells$MetacellsMatrix[ClonePlot,
    ]
    MetaCell_df$CT_bin <- CThripsosObject$Metacells$CT_MetacellsBins[ClonePlot,
    ]
    MetaCell_df$Coordinates <- 1:length(MetaCell_df$CNV)
    MetaCell_df$CT_bin[1] <- 0
    MetaCell_df$CT_bin[2] <- 0.1
    MetaCell_df_filtered <- MetaCell_df
    MetaCell_df_filtered$CT_bin[1] <- 0
    MetaCell_df_filtered$CT_bin[2] <- 0.1
    MetaCell_df_filtered$CNV <- as.numeric(MetaCell_df_filtered$CNV)
    MetaCell_df_filtered$CNV[MetaCell_df_filtered$CNV >=
                               min_cnv_changes] <- min_cnv_changes

    if(score_binary==T)
    {
    MetaCell_df_filtered$CT_bin[MetaCell_df_filtered$CT_bin>0]<-1
    }

    PlotClone <- ggplot(MetaCell_df_filtered, aes(Coordinates,
                                                  CNV)) + geom_point(aes(colour = (CT_bin))) + scale_colour_gradient2(low = "yellow",
                                                                                                                      high = "black", mid = "red", midpoint = (0.5)) +
      ylim(0, 10) + theme_bw() + theme(panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),
                                       axis.text.x = element_blank(), panel.spacing.x = unit(0,
                                                                                             "lines"), panel.border = element_rect(linetype = 3)) +
      facet_grid(. ~ reorder(chrXY, Coordinates), scales = "free",
                 space = "free") + ggExtra::removeGrid() + scale_x_continuous(expand = c(0.01,
                                                                                         0.01))
    plots[[p_i]] <- PlotClone
    p_i = p_i + 1
  }
  AllClones <- gridExtra::grid.arrange(grobs = plots, nrow = rows)
  plot(AllClones)
}


CT_Regions_Metacells<-function(CThripsosObject, Metacell)
{
  #Metacell refers to the metacell index in the matrix e.g 1,2,3..,n
  normalised_ct<-c()
  clone_limits<-c()
  regions<-c()
  
  for(chr_ct in unique(CThripsosObject$Annotations$segment_chromosomes))
  {
    CT_chr<-CThripsosObject$Metacells$CT_MetacellsBins[Metacell,which(CThripsosObject$Annotations$segment_chromosomes==chr_ct)]
    normalised_ct<-c(normalised_ct, CT_chr)
  }
  
  #we normalise the profile to 0,1
  
  # we set a threshold to consider CT+
  CT_thr<-0
  normalised_ct[normalised_ct>CT_thr]<-1
  
  # we take Nans out
  Nans<-which(!is.na(normalised_ct))
  normalised_ct<-normalised_ct[Nans]
  segment_chromosomes_ct<-CThripsosObject$Annotations$segment_chromosomes[Nans]
  segment_coords_ct<-CThripsosObject$Annotations$segment_coords[Nans,]
  
  # we check the profile per chr
  segments_chr_passed<-0
  for(chr_ct in unique(segment_chromosomes_ct))
  {
    # print(paste("chr", chr_ct))
    segments_chr_ct<-which(segment_chromosomes_ct==chr_ct)
    normalised_ct_chr<-normalised_ct[segments_chr_ct]
    segment_chromosomes_ct_chr<-segment_chromosomes_ct[segments_chr_ct]
    segment_coords_ct_chr<-segment_coords_ct[segments_chr_ct,]
    
    compressed<-rle(normalised_ct_chr)
    # compressed$values
    # compressed$lengths
    
    segments_from<-1
    for(i in 1:length(compressed$values))
    {
      if(compressed$values[i]==1)
      {
        # print(i)
        segments_to<-segments_from+compressed$lengths[i]-1
        print(paste("chr", segment_chromosomes_ct_chr[segments_from], segment_coords_ct_chr[segments_from,1], segment_coords_ct_chr[segments_to,2], "segments: ", segments_from, segments_to))
        regions<-rbind(regions, c(segment_chromosomes_ct_chr[segments_from], segment_coords_ct_chr[segments_from,1], segment_coords_ct_chr[segments_to,2]))
        clone_limits<-c(segment_coords_ct_chr[segments_from,1], segment_coords_ct_chr[segments_to,2])
      }
      segments_from<-segments_from+compressed$lengths[i]
    }
    
    segments_chr_passed<-segments_chr_passed + length(segments_chr_ct)
    # plot(normalised_ct_chr, col=as.numeric(segment_chromosomes_ct_chr))
    
  }
  return(regions)
}
