# dyn.load("~/fastPosition.so")

.firstCross <- function(Data, threshold=0, relationship, start=1) {
#Copyright https://github.com/joshuaulrich, under GPLv3
  rel <- switch(relationship[1],
                '>'    =  ,
                'gt'   = 1,
                '<'    =  ,
                'lt'   = 2,
                '=='   =  ,
                'eq'   = 3,
                '>='   =  ,
                'gte'  =  ,
                'gteq' =  ,
                'ge'   = 4,
                '<='   =  ,
                'lte'  =  ,
                'lteq' =  ,
                'le'   = 5)
  .Call('firstCross', Data, threshold, rel, start)
}



fasterCTScore <- function(cell, window_length, min_cnv_changes, min_consec_cnvs, CThripsosObject){

  segment_coords <- apply(CThripsosObject$Annotations$segment_coords,2,as.numeric)

  bin_size = diff(t(segment_coords)[,1])
  if(!all(diff(t(segment_coords))==bin_size)){
    stop("The fast calculation makes an assumption of equal bin sizes.")
  }

  chromatriptic_cell<-FALSE
  chromosomes<- unique(CThripsosObject$Annotations$segment_chromosomes)

  chromatriptic_chromosomes<-c()
  bins_windows<-rep(0, nrow(CThripsosObject$CNVMatrix))
  bins_cnv_changes<-rep(0L, nrow(CThripsosObject$CNVMatrix))


  # chromosome_i <- 2
  for (chromosome_i in chromosomes) {
    # if(chromosome_i == 21){
    #   browser()
    # }
    chromatriptic_chromosome <-0
    chr_windows<-0
    segments_in_chromosome_i=which(CThripsosObject$Annotations$segment_chromosomes==chromosome_i)

    # Window size could be bigger than a chromosome, so we set window size to teh chromosome size.
    if(dplyr::last(segment_coords[segments_in_chromosome_i,2])-dplyr::first(segment_coords[segments_in_chromosome_i,1])<=window_length){
      warning(paste0("Window size is larger than chromosome: ", chromosome_i, "\n"))
      window_length_cur <- dplyr::last(segment_coords[segments_in_chromosome_i,2])-dplyr::first(segment_coords[segments_in_chromosome_i,1])
      # chromatriptic_chromosomes=c(chromatriptic_chromosomes, chromatriptic_chromosome)
      # next
    } else {
      window_length_cur <- window_length
    }

    window_bin_length <- window_length_cur/bin_size

    if(floor(window_bin_length)!=window_bin_length){
      ## TODO: figure out why Gonzalo's implementation is off by one for window length?
      warning("Window is not a multiple of bin size in length. Rounding the window down to the nearest integral bin size.")
      window_bin_length <- floor(window_bin_length)
    }

    window_vector <- rep(1, window_bin_length)



    compressed_cell<-rle(cell[segments_in_chromosome_i])

    windows_to_remove <- which(compressed_cell$lengths<=min_consec_cnvs)
    if(length(windows_to_remove)==length(compressed_cell$lengths)){
      stop(paste0("Error: Chromosome:", chromosome_i, "has no segments longer than min_consec_cnvs. Check that the parameter value is sane, and that CNVs segmentation has no artefacts. Also check that your bins are sorted numerically, and not alphabetically."))
    }

    for(windw in windows_to_remove){
      if(windw==1){
        ## Special case, we need to find the first valid segment in the chromosome and set it to that value, so the rest of the segments in front will be also correct.
        compressed_cell$values[windw] = compressed_cell$values[which(compressed_cell$lengths>=min_consec_cnvs)[[1]]]
      } else {
        compressed_cell$values[windw] = compressed_cell$values[windw-1]
      }
    }

    ## double check that this works properly for non-even segment/bin sizes
    cell[segments_in_chromosome_i] <- inverse.rle(compressed_cell)

    compressed_cell <- rle(inverse.rle(compressed_cell))

    if (length(compressed_cell$values)<min_cnv_changes) {
      # print("skipping chr...")
      chromatriptic_chromosomes=c(chromatriptic_chromosomes, chromatriptic_chromosome)
      next
    }

    breakpoint_idx <- cumsum(compressed_cell$lengths)[-length(compressed_cell$lengths)] ## Need to remove the last value as there is no breakpoint there.

    segment_coords_chromosome = segment_coords[segments_in_chromosome_i,]

    ## This is the right/end of the consective same copy number area, but therefore left of the breakpoint in question (since a breakpoint is in between genomic loci)
    breakpoint_locs_left <- as.numeric(segment_coords_chromosome[breakpoint_idx,2])

    breakpoints <- seq_along(breakpoint_locs_left)


    ## R doesn't allow you to jump ahead in for loops, so this is a while loop.
    i = 1
    while(i <= (length(breakpoint_locs_left) - min_cnv_changes + 1)){
      # print(paste0("i at beginning of loop: ",i))
      ## If the next condition is true, we have a chromothripsis event.
      if(diff(breakpoint_locs_left[c(i,i+min_cnv_changes-1)])<=window_length_cur){ # include the right edge of last breakpoint into comparison with the length of my window

        ## Now, we need to find the first window in which this chromothripsis event occurs. We do this "generously", that is, we require any section of the bin
        ## to intersect the window to include it in the event. This can trivially be made conservative by instead requiring the full bins to be inside the
        ## window.

        active_breakpoints <- breakpoints[seq(i, i+min_cnv_changes - 1)]


        ## sequence generation in R can be frustrating, so this is a while loop.
        j <- i + 1
        while(j <= length(breakpoint_locs_left) - min_cnv_changes + 1){
          if(diff(breakpoint_locs_left[c(j,j+min_cnv_changes-1)])<=window_length_cur){
            active_breakpoints <- c(active_breakpoints, breakpoints[j+min_cnv_changes-1])
          }
          j <- j + 1
        }


        active_breakpoints_locleft <- breakpoint_locs_left[active_breakpoints]

        distance_between_ends <- last(active_breakpoints_locleft) - dplyr::first(active_breakpoints_locleft)
        bins_between_ends <- breakpoint_idx[last(active_breakpoints)] - breakpoint_idx[dplyr::first(active_breakpoints)] - 1

        leftmost_location_with_ct <- max(active_breakpoints_locleft[min_cnv_changes] - window_length_cur + 1, 1)
        rightmost_location_with_ct <- min(active_breakpoints_locleft[length(active_breakpoints_locleft) - min_cnv_changes + 1] + window_length_cur - 1, as.numeric(tail(segment_coords_chromosome[,1],1)))

        leftmost_bin_with_ct <- Position(function(y) y >= leftmost_location_with_ct, as.numeric(segment_coords_chromosome[,2])) # this is the "anti-conservative" decision
        rightmost_bin_with_ct <- Position(function(y) y <= rightmost_location_with_ct, as.numeric(segment_coords_chromosome[,1]), right = TRUE)


        ## TODO:: Experiment if this can be made even faster by not doing the convolution (or is FFT fast enough that this doesn't really matter?)
        # ## taking care of endpoints being defined by chromosome end. These endpoints are included in the repetitive region
        # left_end_of_repeating_region <- leftmost_bin_with_ct + window_bin_length
        # right_end_of_repeating_region <- Position(function(y) y <= leftmost_location_with_ct + window_length, as.numeric(segment_coords_chromosome[,2]), right=TRUE)
        #

        # num_bins_before_repeating_region <- left_end_of_repeating_region - leftmost_bin_with_ct  # + 1 because the leftmost bin is included
        # num_bins_after_repeating_region <- rightmost_bin_with_ct - right_end_of_repeating_region # the location of the last active bin is still in between the breakpoints
        #
        # number_in_repeating_region <- min(num_bins_before_repeating_region,num_bins_after_repeating_region) + 1

        # bins_windows[segments_in_chromosome_i][seq(leftmost_bin_with_ct, left_end_of_repeating_region-1)] <-
        #   bins_windows[segments_in_chromosome_i][seq(leftmost_bin_with_ct, left_end_of_repeating_region-1)] + seq(1,num_bins_before_repeating_region)
        #
        # bins_windows[segments_in_chromosome_i][seq(left_end_of_repeating_region, right_end_of_repeating_region)] <-
        #   bins_windows[segments_in_chromosome_i][seq(left_end_of_repeating_region, right_end_of_repeating_region)] + number_in_repeating_region
        #
        # bins_windows[segments_in_chromosome_i][seq(right_end_of_repeating_region+1,rightmost_bin_with_ct)] <-
        #   bins_windows[segments_in_chromosome_i][seq(right_end_of_repeating_region+1,rightmost_bin_with_ct)] + seq(num_bins_after_repeating_region,1)
        # chromatriptic_chromosome <- chromatriptic_chromosome + number_in_repeating_region + 1

        ## TODO: This convolution is now the slowest part. I think I can simplify this even further, as a general convolution is not actually necessary here, I know these filters will always be all 1s,
        ## I should be able to compute the result from the length of the filters. However, I think this is efficient enough for now.

        conv_res = convolve(rep(1,rightmost_bin_with_ct - window_bin_length - leftmost_bin_with_ct + 2), window_vector, type = "o") ## + 2 because we want to include the end bins, as well as only a single fitting CT window is still ok.
        chromatriptic_chromosome <- chromatriptic_chromosome + length(conv_res) - window_bin_length + 1
        bins_windows[segments_in_chromosome_i][seq(leftmost_bin_with_ct, rightmost_bin_with_ct)] <- bins_windows[segments_in_chromosome_i][seq(leftmost_bin_with_ct, rightmost_bin_with_ct)] + conv_res
        ## Skip to the end of the currently active breakpoints
        i <-  last(active_breakpoints)
        # print(paste0("i at end of loop: ", i))
      }
      i <- i + 1
    }

    ## bin cnv changes is done on all breakpoints:
    ## TODO:: this does not currently agree with the reference implementation, but I am not convinced the reference is doing the right thing, because I think this does what I expect.
    ii <- dplyr::first(breakpoints)

    while(ii <= last(breakpoints)){
      # jj <- ii + min_cnv_changes

      ## Finding the rightmost bin that breakpoint ii has influence on.
      rightmost_location_with_ct_current <- min(breakpoint_locs_left[ii] + window_length_cur - 1, as.numeric(tail(segment_coords_chromosome[,2],1))-1)   ## We want to fully incorporate the breakpoint here
      rightmost_bin_with_ct_current <- .firstCross(segment_coords_chromosome[,2], threshold = rightmost_location_with_ct_current, relationship = "gt")

      jj <- ii+1
      while(jj <= last(breakpoints)){
        ## are ii and jj close enough together?
        if(diff(breakpoint_locs_left[c(ii,jj)]) <= window_length_cur){
          leftmost_location_with_ct_current <- max(breakpoint_locs_left[jj] - window_length_cur + 1, 1)  ## We want to fully incorporate the breakpoint here
          leftmost_bin_with_ct_current <- .firstCross(segment_coords_chromosome[,2], threshold = leftmost_location_with_ct_current, relationship = "gteq")

          bins_cnv_changes[segments_in_chromosome_i][leftmost_bin_with_ct_current:rightmost_bin_with_ct_current] <- pmax(jj - ii + 1,bins_cnv_changes[segments_in_chromosome_i][leftmost_bin_with_ct_current:rightmost_bin_with_ct_current])
        }
        jj <- jj + 1
      }
      ii <- ii + 1
    }


    # chr_windows <- last_segment_in_valid_window <- Position(function(x) tail(as.numeric(segment_coords_chromosome[,2]),1) - window_length > x, as.numeric(segment_coords_chromosome[,1]), right=TRUE)
    chr_windows <- length(segments_in_chromosome_i) - window_bin_length + 1
    if(chromatriptic_chromosome >0) {
      # print(paste0("Chromosome: ", chromosome_i, " Number of CT windows: ", chromatriptic_chromosome, " Total number of windows: ", chr_windows))
      chromatriptic_chromosome <- chromatriptic_chromosome/chr_windows
      bins_windows[segments_in_chromosome_i]<-bins_windows[segments_in_chromosome_i]/chr_windows
    }
    chromatriptic_chromosomes=c(chromatriptic_chromosomes, chromatriptic_chromosome)

  }

  # print(chromatriptic_chromosomes)
  return(list(chromatriptic_chromosomes=chromatriptic_chromosomes, bins_windows=bins_windows, bins_cnv_changes=bins_cnv_changes))
}
