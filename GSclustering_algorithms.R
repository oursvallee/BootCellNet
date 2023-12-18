# gene-set-based refinement of cell clustering
#
# see the followings
# C:/work/20230212_cellstrap/20231001_scrna_test_8_covid19/NatMed_2020_Liao/20231016_analysis_Liao-CoV_1.R
# C:/work/20230212_cellstrap/20231001_scrna_test_8_covid19/NatMed_2020_Liao/20231016_analysis_Liao-CoV_2.R
# C:/work/20230212_cellstrap/20231001_scrna_test_8_covid19/NatMed_2020_Liao/20231016_analysis_Liao-CoV_3.R
#
# for usage, see C:/work/20230212_cellstrap/20231001_scrna_test_8_covid19/PBMC/20231115_pbmc_analysis_1.R
#
# 11/16/2023; 12/8/2023
########################################################

#=======================================================
# averaging the data based on the given gene set
# gset: a set of genes to be checked
# sobj: Seurat object; active.ident should be set before doing this

avgGS <- function(gset, sobj){
  
  ident <- sort(unique(sobj@active.ident))
  data <- (sobj[[sobj@active.assay]]@data)[gset,]
  
  data.avg <- matrix(nrow=length(gset), ncol=length(ident))
  
  # for(i in sort(unique(so@active.ident))){print(sum(so@active.ident==i))}
  # for(i in sort(unique(so@active.ident))){print(c(i,sum(so@active.ident==i)))}
  
  for(i in 1:length(ident)){
    data.curr <- data[,sobj@active.ident==ident[i]]
    if(sum(sobj@active.ident==ident[i])>1){
      data.avg[,i] <- rowMeans(data.curr)  
    }
    else{
      data.avg[,i] <- data.curr
    }
  }
  
  rownames(data.avg) <- gset
  colnames(data.avg) <- ident
  
  return(data.avg)
}



#=======================================================
# put new clusters onto the Seurat object 
# cf 20231124_pbmc_analysis_3.R
# sobj: Seurat object; active.ident should be set before doing this
# clst: clusters to be put
# 12/12/2023: unfactorize the vector
# make it faster

# wrapper
newclst <- function(sobj, clst, clst.name=NULL){
  
  clst.new <- map.newclst(sobj, clst)
  levels <- unique(clst)
  
  clst.new <- factor(clst.new, levels=levels)
  if(is.null(clst.name)){
    clst.name <- "MDS_clst"
  }
  sobj[[clst.name]] <- clst.new
  
  return(sobj)
}

#----------------------------
# auxiliary function to map the clusters
map.newclst <- function(sobj, clst) {
  active.ident <- sobj@active.ident
  unique.active.ident <- unique(active.ident)
  
  # Use match and vectorized operations for faster performance
  clst.new <- clst[match(active.ident, unique.active.ident)]
  names(clst.new) <- names(active.ident)
  
  return(clst.new)
}


#=======================================================
# calculate adjusted Rand index
# give two clusters, and return the index

# wrapper
AdjRandInx <- function(sobj, clst1, clst2){
  clst1.so <- map.newclst(sobj, clst1)
  clst2.so <- map.newclst(sobj, clst2)
  
  return(calc.ari(clst1.so, clst2.so))
}

# auxiliary func
calc.ari <- function(x, y){
  
  cont.tab <- table(x, y)
  cont.tab.bco <- sum(cont.tab*(cont.tab-1)/2)
  
  rsum <- rowSums(cont.tab); csum <- colSums(cont.tab)
  rsum.bco <- sum(rsum*(rsum-1))/2; csum.bco <- sum(csum*(csum-1))/2
  tot.bco <- sum(rsum)*(sum(rsum)-1)/2
  
  ari <- (cont.tab.bco-rsum.bco*csum.bco/tot.bco)/((rsum.bco+csum.bco)/2-rsum.bco*csum.bco/tot.bco)
  
  return(ari)
  
}

#=======================================================
# resampling/shuffling to calculate significance for ARI
# sobj: Seurat object, where clustering 1/2 was done
# clst1: fixed
# clst2: shuffled
# N: number of resampling

rs.ARI <- function(sobj, clst1, clst2, N=1000){
  
  clst.fixed <- map.newclst(sobj, clst1)
  unique.sobj.clst <- as.vector(sobj@active.ident)
  
  aris <- rep(0, N)
  
  for(j in 1:N){
    clst2.shuffled <- clst2[sample(1:length(clst2))]
    names(clst2.shuffled) <- names(clst2)

    clst.shuffled <- map.newclst(sobj, clst2.shuffled)
    aris[j] <- calc.ari(clst.fixed, clst.shuffled)
  }
  
  return(aris)
}




###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# map.newclst <- function(sobj, clst){
#   
#   active.ident <- sobj@active.ident
#   clst.new <- as.vector(active.ident)
#   # names(clst.new) <- names(active.ident)
#   # n.clst <- max(clst)
#   
#   
#   unique.active.ident <- sort(unique(active.ident))
#   for(i in 1:length(unique.active.ident)){
#     clst.new[active.ident==unique.active.ident[i]] <- clst[i]
#   }
#   
#   return(clst.new)
#   
# }