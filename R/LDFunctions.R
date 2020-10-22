# Split output of VCFtools, or similar program,
# into chromsome-by-chromosome list
# INPUT
# LD = matrix from vcftools, with character column "CHR", and numeric colums "POS1", "POS2", "N_IND", "R.2"
# nchrom = number of chromosomes represented in the LD file
# sep = character to split the CHR values at to get chromosome names
splitByChrom <- function(LD, nchrom = 9, sep = "_") {
  LD.by.chroms <- vector(mode = "list", length = nchrom)
  
  if(sep != ""){
    allmarkers <- sapply(strsplit(LD$CHR, sep), "[[", 2)
    chroms <- unique(allmarkers)
  }  else{
    allmarkers <- LD$CHR
    chroms <- unique(allmarkers)
  }
  # Split into chromsome-by-chromosome list
  for(i in 1:nchrom){
    chrom <- chroms[i]
    pos <- which(allmarkers == chrom)
    LD.by.chroms[[i]] <- LD[pos,]
  }
  return(LD.by.chroms)
}



# Get the distance between pairs of markers of VCFtools, or similar program,
# into chromsome-by-chromosome list
# INPUT
# LD = list of dataframes, one for each chromsome
# Each element of the list should have a character column "CHR", and numeric colums "POS1", "POS2", "N_IND", "R.2"
calcDist <- function(LD) {
  
  library(tictoc)
  data <- vector(mode = "list", length = length(LD)+1)
  
  for(i in 1:length(LD)){
    
    print(paste("Processing pairs of markers on chromosome:", i))
    r2.file <- LD[[i]]
    ans <- matrix(nrow = nrow(r2.file), ncol = 2)
    colnames(ans) <- c("dist", "r2")
    
    print(paste(nrow(r2.file), "marker pairs to process"))
    tic()
    for(j in 1:nrow(r2.file)){
      
      if(j %% 100000 == 0){
        print(paste("Processing marker pair", j, "of", nrow(r2.file)))
      }
      pair <- r2.file[j,]
      dist <- abs(pair$POS1 - pair$POS2)
      r2 <- pair$R.2
      ans[j,1] <- dist
      ans[j,2] <- r2
    }
    toc()
    
    data[[i]] <- ans
  }
  
  # Add back the full dataset and save
  data[[length(LD)+1]] <- do.call("rbind", data)
  listNames <- as.character(seq(1:length(data)))
  listNames[length(LD)+1] <- "all"
  names(data) <- listNames
  
  return(data)
}



# Get the mean of each bin of SNPs, with size determined by vcftool call, and specified by windowSize
# INPUT
# LD = list of dataframes, one for each chromsome
# Each element of the list should have a character column "CHR", and numeric colums "POS1", "POS2", "N_IND", "R.2"
slidingWindow <- function(LD, windowSize = 100) {
  library(tictoc)
  data <- vector(mode = "list", length = length(LD)+1)
  
  for(i in 1:length(LD)){
    tic()
    
    r2.file <- LD[[i]]
    data.sub <- matrix(data = NA, nrow = ceiling(nrow(r2.file)/windowSize), ncol = 3)
    colnames(data.sub) <- c("POS", "R2", "DIST")
    
    bins <- seq(from=1,to=nrow(r2.file),by=windowSize)
    
    print(paste("Total number of windows on chromsome", i, ":", length(bins)))
    
    for(j in 1:(length(bins)-1)){
      r2.sub <- r2.file[c(bins[j]:(bins[j]+(windowSize - 1))),]
      r2.sub$DIST <- r2.sub$POS2 - r2.sub$POS1
      data.sub[j, 1] <- mean(r2.sub$POS1)
      data.sub[j, 2] <- mean(r2.sub$R.2)
      data.sub[j, 3] <- mean(r2.sub$DIST)
    }
    
    # Last bin
    j <- j+1
    r2.sub <- r2.file[c(bins[j]:nrow(r2.file)),]
    r2.sub$DIST <- r2.sub$POS2 - r2.sub$POS1
    data.sub[j, 1] <- mean(r2.sub$POS1)
    data.sub[j, 2] <- mean(r2.sub$R.2)
    data.sub[j, 3] <- mean(r2.sub$DIST)
    
    data.sub <- as.data.frame(data.sub)
    
    data.sub$chrom <- i
    data.sub <- data.sub[,c(4,1,2,3)]
    colnames(data.sub) <- c("chr", "bp", "r2", "dist")
    data[[i]] <- data.sub
    toc()
  }
  
  
  data[[length(LD)+1]] <- do.call("rbind", data)
  listNames <- as.character(seq(1:length(data)))
  listNames[length(LD)+1] <- "all"
  names(data) <- listNames
  return(data)
}






# As above, but with non-overlapping windows -- i.e., the start point of the next window
# is the position of the SNP farthest away from the previous SNP's 100-SNP bin.  I think this is what 
# Shelby asked me to compare.
nonOverlappingWindow <- function(LD, windowSize = 100) {
  library(tictoc)
  data <- vector(mode = "list", length = length(LD)+1)
  
  for(i in 1:length(LD)){
    tic()
    
    r2.file <- LD[[i]]
    
    data.sub <- matrix(data = NA, nrow = ceiling(nrow(r2.file)/windowSize), ncol = 2)
    colnames(data.sub) <- c("POS", "R2")
    
    print(paste("Processing chromsome", i))
    
    endPoint <- 0
    start <- which(r2.file$POS1 > endPoint)[1]
    end <- start+(windowSize - 1)
    
    j <- 1
    
    while(is.na(r2.file[end, 3]) == F)  {
      r2.sub <- r2.file[c(start:end),]
      
      data.sub[j, 1] <- mean(r2.sub$POS1)
      data.sub[j, 2] <- mean(r2.sub$R.2)
      j <- j + 1
      
      endPoint <- r2.sub[windowSize,3]
      start <- which(r2.file$POS1 > endPoint)[1]
      end <- start+(windowSize - 1)
    }
    data.sub <- as.data.frame(data.sub)
    data.sub <- data.sub[is.na(data.sub$POS) == F,]
    
    data.sub$chrom <- i
    data.sub <- data.sub[,c(3,1,2)]
    colnames(data.sub) <- c("chr", "bp", "r2")
    data[[i]] <- data.sub
    toc()
  }
  
  # Add in the entire dataset
  data[[length(LD)+1]] <- do.call("rbind", data)
  
  # Name the list elements after the chroms
  listNames <- as.character(seq(1:length(data)))
  listNames[length(LD)+1] <- "all"
  names(data) <- listNames
  
  return(data)
}






# Fit exponential decay to the LD pair data
# into chromsome-by-chromosome list
# INPUT
# LD = list of dataframes, one for each chromsome
# Each element of the list should have a numeric column "dist" and "r2"
fitDecay <- function(LD, chrom = 10) {
  markerPairs <- as.data.frame(LD[[chrom]])
  
  powers <- seq(from=0,to=6,by=0.10)
  markerPairs$group <- cut(markerPairs$dist, breaks=10^powers)
  markerPairs <- markerPairs[is.na(markerPairs$group)==F,]
  markerPairs <- tapply(markerPairs$r2, markerPairs$group, mean, na.rm = TRUE)
  
  markerPairs <- data.frame(group = names(markerPairs), r2 = markerPairs)
  markerPairs <- markerPairs %>% mutate(start=as.integer(str_extract(str_replace_all(group,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                            end=as.integer(str_extract(str_replace_all(group,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                            mid=start+((end-start)/2))
  markerPairs <- markerPairs[,c(5,2)]
  fit_1M <- nls(r2 ~ SSasymp(log10(mid), yf, y0, log_alpha), data = markerPairs)
  
  summaryTable <- summary(fit_1M)
  params <- summaryTable$coefficients[,1]
  
  Asym <- params[1]
  R0 <- params[2]
  lrc <- params[3]
  
  x <- ((Asym - R0)/(Asym - 0.1))^(exp(-lrc)*log(10))
  x0.2 <- ((Asym - R0)/(Asym - 0.2))^(exp(-lrc)*log(10))
  
  p_1M_log <- qplot(mid, r2, data = augment(fit_1M), size=I(2.5)) +
    # geom_point(data = data[[10]], aes(x=dist, y=r2), size=2, shape=23)+
    geom_line(aes(y = .fitted),size=1.5,alpha=1,colour="red") +
    geom_segment(aes(x = x, xend = x, y = 0, yend = 0.1), color = "green", size = 1.5) +
    geom_segment(aes(x = x0.2, xend = x0.2, y = 0, yend = 0.2), color = "blue", size = 1.5) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color="blue", size = 1.5) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color="green", size = 1.5) +
    labs(x="Distance (bp)",y=expression("Mean LD"~(r^{2}))) +
    scale_x_continuous(trans="log10", breaks = c(1*10^0,1*10^1,1*10^2,x0.2,1*10^3,1*10^4,x,1*10^5,1*10^6), labels=c("1","10","100","796", "1k","10k","19.7k", "100k","1M"))+#breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6,10*10^6,12*10^6),labels=c("0","2","4","6","8","10","12"))+
    scale_y_continuous(expand = c(0,0), limits = c(0,0.7))+#breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6,10*10^6,12*10^6),labels=c("0","2","4","6","8","10","12"))+
    theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ans <- list(
    modelFit <- fit_1M,
    summaryTable <- summaryTable,
    r2_0.10 <- x,
    plot <- p_1M_log)
  
  names(ans) <- c("nlsModel", "nlsSummary", "r2cutoff", "ggplot")
  return(ans)
}




## All credit to user "rmf": https://www.biostars.org/p/347796/
plotPairwiseLD <- function(dfr, chr, xlim=c(NA,NA), ylim=c(NA,NA), minr2) {
  if(missing(dfr)) stop("Input data.frame 'dfr' missing.")
  ld <- dfr
  if(!missing(chr)) {
    ld <- filter(ld,CHROM_A==get("chr") & CHROM_B==get("chr"))
  }
  ld <- filter(ld,POS_A<POS_B)
  
  if(!missing(minr2)) {
    ld <- filter(ld,R2>get("minr2"))
  }
  
  ld <- ld %>% arrange(R2)
  
  ld$x <- ld$POS_A+((ld$POS_B-ld$POS_A)/2)
  ld$y <- ld$POS_B-ld$POS_A
  ld$r2c <- cut(ld$R2,breaks=seq(0,1,0.2),labels=c("0-0 - 0.2","0.2 - 0.4",
                                                   "0.4 - 0.6","0.6 - 0.8",
                                                   "0.8 - 1.0"))
  ld$r2c <- factor(ld$r2c,levels=rev(c("0-0 - 0.2","0.2 - 0.4",
                                       "0.4 - 0.6","0.6 - 0.8",
                                       "0.8 - 1.0")))
  
  ggplot(ld,aes(x=x,y=y,col=r2c))+
    geom_point(shape=20,size=0.1,alpha=0.8)+
    scale_color_manual(values=c("#ca0020","#f4a582","#d1e5f0","#67a9cf","#2166ac"))+
    scale_x_continuous(limits=xlim)+
    scale_y_continuous(limits=ylim)+
    guides(colour=guide_legend(title="R2",override.aes=list(shape=20,size=8)))+
    labs(x="Chromosome 8 (bp)",y="")+
    theme_bw(base_size=14)+
    theme(panel.border=element_blank(),
          axis.ticks=element_blank()) %>%
    return()
}



