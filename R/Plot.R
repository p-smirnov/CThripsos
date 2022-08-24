plot_MetacellsCT<-function (CThripsosObject, score_binary=T, rows=NULL, plotvar="CT", max_cnv=CThripsosObject$Metacells$min_cnv_changes)
{
  if( is.null(CThripsosObject$Metacells$CT_MetacellsBins))
  {
    print("No CT data was found. Did you run Calculate_CT_Metacells()?")
  }
  else{


  # max_cnv is the upper limit to the color scale
  library(ggplot2)

  if(is.null(rows))
  {
    rows<-nrow(CThripsosObject$Metacells$MetacellsMatrix)
  }

  MetaCell_All_df <- c()
  plots <- list()
  p_i = 1
  for (ClonePlot in unique(CThripsosObject$Metacells$MetacellsClusters[,2]))
  {
    MetaCell_df <- as.data.frame(cbind(rep(ClonePlot, length(CThripsosObject$Annotations$segment_chromosomes)),
                                       colnames(CThripsosObject$Metacells$MetacellsMatrix)))

    colnames(MetaCell_df) <- c("Clone", "Chromosome")
    MetaCell_df$chrXY <- unlist(lapply(strsplit(MetaCell_df$Chromosome,
                                                ":"), "[[", 1))
    MetaCell_df$start <- unlist(lapply((strsplit(unlist(lapply(strsplit(MetaCell_df$Chromosome, "-"), "[[", 1)), ":")), "[[", 2))
    MetaCell_df$end <- unlist(lapply(strsplit(MetaCell_df$Chromosome, "-"), "[[", 2))
    MetaCell_df$CNV <- CThripsosObject$Metacells$MetacellsMatrix[ClonePlot, ]
    MetaCell_df$CT_bin <- CThripsosObject$Metacells$CT_MetacellsBins[ClonePlot,]
    MetaCell_df$CT_bin_cnvmax <- CThripsosObject$Metacells$CT_MetacellsBinsMaxCN[ClonePlot,]
    MetaCell_df$CT_bin_cnvmax[which(MetaCell_df$CT_bin_cnvmax >  max_cnv)] <- max_cnv
    MetaCell_df$CT_bin_cnvmax[1] <- 0
    MetaCell_df$CT_bin_cnvmax[2] <- min_cnv_changes

    MetaCell_df$Coordinates <- 1:length(MetaCell_df$CNV)
    MetaCell_df$CT_bin[1] <- 0
    MetaCell_df$CT_bin[2] <- 0.1
    MetaCell_df_filtered <- MetaCell_df
    MetaCell_df_filtered$CT_bin[1] <- 0
    MetaCell_df_filtered$CT_bin[2] <- 0.1
    MetaCell_df_filtered$CNV <- as.numeric(MetaCell_df_filtered$CNV)
    MetaCell_df_filtered$CNV[MetaCell_df_filtered$CNV >= min_cnv_changes] <- min_cnv_changes

    if(score_binary==T)
    {
      MetaCell_df_filtered$CT_bin[MetaCell_df_filtered$CT_bin>0]<-1
    }

    if(plotvar=="CT")
    {
      PlotClone <- ggplot(MetaCell_df_filtered, aes(Coordinates, CNV)) + ggrastr::rasterise(geom_point(aes(colour = (CT_bin))), dpi=300) +
        scale_colour_gradient2(low = "orange", high = "black", mid = "red", midpoint = (0.5)) +
        ylim(0, 10) + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x = unit(0, "lines"), panel.border = element_rect(linetype = 3)) +
        facet_grid(. ~ reorder(chrXY, Coordinates), scales = "free", space = "free") + ggExtra::removeGrid() + scale_x_continuous(expand = c(0.01, 0.01))
    }else if(plotvar=="MaxCNV"){
      PlotClone <- ggplot(MetaCell_df_filtered, aes(Coordinates, CNV)) + ggrastr::rasterise(geom_point(aes(colour = (CT_bin_cnvmax))), dpi=300) +
        scale_colour_gradient2(low = "orange", high = "black", mid = "red", midpoint = (5)) +
        ylim(0, 10) + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x = unit(0, "lines"), panel.border = element_rect(linetype = 3)) +
        facet_grid(. ~ reorder(chrXY, Coordinates), scales = "free", space = "free") + ggExtra::removeGrid() + scale_x_continuous(expand = c(0.01, 0.01))
    }

    plots[[p_i]] <- PlotClone
    p_i = p_i + 1
  }
  AllClones <- gridExtra::grid.arrange(grobs = plots, nrow = rows)

  return(AllClones)
  }
}

