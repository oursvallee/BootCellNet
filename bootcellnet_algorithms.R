# BootCellNet algorithms
# 12/6/2023
#
# summary of 
# C:/work/20230212_cellstrap/20230620_scrna_test_5/20230620_1_test.R
#
# for usage, see 
# C:/work/20230212_cellstrap/20230620_scrna_test_5/20230620_1_test.R
# C:/work/20230212_cellstrap/20231001_scrna_test_8_covid19/PBMC/20231108_pbmc_analysis_2.R
#
########################################################


##################################################
# procdure
# (1) K-nn smoothing
# (2) "cellstrap": (i) select N cells for M times, (ii) calculate GRN, (iii) write down the data into a file
# (3) aggregate the data; count TFs, count edges
# (4) draw graph: weighted nodes and edges?

#===================================
# (1) K-nn smoothing **possible improvement here**
# cf C:/work/20230212_cellstrap/20230710_scrna_test_6/20230710_1_test.R
# adaptive thresholding w/ Otsu method has been implemented, but does not affect much...
# gzip (7/25/2023 added)
# avoid data retrieval from dgCMatrix, faster than original (10/24/2023)

knnsmo <- function(data, umap, K=50, theta=1, file, progress=T){
  require(Matrix)
  M0 <- dim(data)[1]
  N0 <- dim(data)[2] # cell number
  
  # data.p <- data@p
  # data.x <- data@x
  # data.i <- data@i
  
  x.cs <- integer(length=N0*M0)
  i.cs <- integer(length=N0*M0)
  j.cs <- integer(length=N0*M0)
  nonz <- 0
  nonz.prev <- 1
  
  for(j in 1:N0){
    # print(j)
    umap.curr <- umap[j,]
    
    dists <- apply(umap, 1, function(x) sum((x-umap.curr)^2))
    dists[j] <- Inf # remove the origin
    
    cs.set <- vector(length = K)
    # data.set <- matrix(nrow=M0, ncol=K)
    cs.vec <- vector(length = M0)
    
    for(k in 1:K){
      k.min <- which.min(dists)
      cs.set[k] <- k.min
      
      # directly access to the sparse data...
      p.k <- data@p[k.min]
      p.k2 <- data@p[k.min+1]
      inx.k <- data@i[(p.k+1):p.k2]+1
      x.k <- data@x[(p.k+1):p.k2]
      
      cs.vec[inx.k] <- cs.vec[inx.k]+x.k 
      dists[k.min] <- Inf
    }
    # cs.set <- order(dists)[2:(K+1)]
    # data.set[is.na(data.set)] <- 0
    cs.vec <- cs.vec/K
    cs.vec[cs.vec<theta] <- 0
    cs.nonz <- cs.vec>0
    n.nonz <- sum(cs.nonz)
    if(n.nonz<1){
      next
    }
    nonz <- nonz+n.nonz
    
    i.cs[nonz.prev:nonz] <- (1:M0)[cs.nonz]
    j.cs[nonz.prev:nonz] <- rep(j,n.nonz)
    x.cs[nonz.prev:nonz] <- floor(cs.vec[cs.nonz])
    
    nonz.prev <- nonz+1
    
    if(progress){
      if(j%%100==0){
        cat(paste0(j, " of ", N0, " done\n"))
      }
    }
    
  }
  
  i.cs <- i.cs[i.cs>0]; i.cs <- c(i.cs, M0)
  j.cs <- j.cs[j.cs>0]; j.cs <- c(j.cs, N0)
  x.cs <- x.cs[x.cs>0]; x.cs <- c(x.cs, 0)
  
  data.cs <- sparseMatrix(i=i.cs, j=j.cs, x=x.cs)
  
  rownames(data.cs) <- rownames(data)
  colnames(data.cs) <- colnames(data.cs)
  
  if(progress){
    cat(paste0("Writing file ", file, "...\n"))
  }
  # writeMM(data.cs, file=file)
  writeMMgz(data.cs, file)
  return(data.cs)
}


# usage
# data.kns <- knnsmo(data=data, umap=umap.coord, file="20230620_test_kns.mtx.gz")

# # you can gzip the mtx file
# data.kns <- readMM("20230620_test_kns.mtx.gz")
# rownames(data.kns) <- rownames(data)
# colnames(data.kns) <- colnames(data)


#===================================
# (1) K-nn smoothing **possible improvement here**
# gzip (7/25/2023 added)
# avoid data retrieval from dgCMatrix, faster than original (10/24/2023)
# writing directly into a gz file (11/1/2023)

knnsmo.dir <- function(data, umap, K=50, theta=1, file, progress=T){
  require(Matrix)
  
  # prepare file
  GZF <- gzfile(file, "at")
  
  M0 <- dim(data)[1]
  N0 <- dim(data)[2] # cell number
  
  # data.p <- data@p
  # data.x <- data@x
  # data.i <- data@i
  
  nonz <- 0
  nonz.prev <- 1
  
  for(j in 1:N0){
    # print(j)
    umap.curr <- umap[j,]
    
    dists <- apply(umap, 1, function(x) sum((x-umap.curr)^2))
    dists[j] <- Inf # remove the origin
    
    cs.set <- vector(length = K)
    cs.vec <- vector(length = M0)
    
    for(k in 1:K){
      k.min <- which.min(dists)
      cs.set[k] <- k.min
      
      # directly access to the sparse data...
      p.k <- data@p[k.min]
      p.k2 <- data@p[k.min+1]
      inx.k <- data@i[(p.k+1):p.k2]+1
      x.k <- data@x[(p.k+1):p.k2]
      
      cs.vec[inx.k] <- cs.vec[inx.k]+x.k 
      dists[k.min] <- Inf
    }
    
    cs.vec <- cs.vec/K
    cs.vec[cs.vec<theta] <- 0
    cs.nonz <- cs.vec>0
    n.nonz <- sum(cs.nonz)
    if(n.nonz<1){
      next
    }
    
    i.cs <- (1:M0)[cs.nonz]
    j.cs <- rep(j,n.nonz)
    x.cs <- floor(cs.vec[cs.nonz])
    
    for(l in 1:n.nonz){
      writeLines(sprintf("%s %s %s", i.cs[l], j.cs[l], x.cs[l]), GZF)
    }
    
    if(progress){
      if(j%%100==0){
        cat(paste0(j, " of ", N0, " done\n"))
      }
    }
    
    nonz <- nonz+n.nonz
    
  }
  
  close(GZF)
  
  headerfile <- paste0("header_",file)
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate real general"),
      sprintf("%s %s %s", M0, N0, nonz)
    ),
    gzfile(headerfile)
  )
  
  filewheader <- paste0("wh_",file)
  if(progress){
    cat("Please cat the header + the matrix market file by running\n")
    cat(paste0("cat ", headerfile, " ", file, " >", filewheader,"\n"))
  }
}

# usage
# knnsmo.dir(data=data, umap=umap.coord, file="20230620_test_kns.mtx.gz")
# cat header_20230620_test_kns.mtx.gz 20230620_test_kns.mtx.gz > wh_20230620_test_kns.mtx.gz

# # you can gzip the mtx file
# data.kns <- readMM("wh_20230620_test_kns.mtx.gz")
# rownames(data.kns) <- rownames(data)
# colnames(data.kns) <- colnames(data)




#===================================
# (2) "BootCellNet": (i) select N cells for M times, (ii) calculate GRN, (iii) write down the data into a file
# goi should be given in a same format as rownames of data
# 11/10/2023: changed the name; added "shuffle" option 
# 11/13/2023: "shuffle" can be either "all", "goi", or "none"
# 11/14/2023: "nested" version, takes more time but can get lots of bootstrap samples!!
# 11/14/2023: "nested" loop is not necessarily run on here; "shuffle" here is all or none.
# 12/6/2023: avoid data retrieval from dgCMatrix, faster than original
# original: with Jorstad-human data and N=200, takes 3 min/specimen...
# new: takes < 10sec with the same, takes 1 min/specimen with the data and N=20000

# data: should be dgCMatrix
# N: number of cells to sample
# M: number of iteration

bootcellnet <- function(data, N=5000, M=100, goi=NULL, 
                        directory, suffix=NULL, 
                        nfeatures=2000, selection.method="vst", 
                        estimator="mi.mm", disc = "equalwidth", 
                        shuffle=FALSE, verbose=FALSE){
  require(Seurat)
  require(minet)
  
  if(!file.exists(directory)){
    dir.create(directory)
  }
  
  G0 <- dim(data)[1] # gene number
  N0 <- dim(data)[2] # cell number
  
  # data.p <- data@p
  # data.x <- data@x
  # data.i <- data@i
  
  for(m in 1:M){
    
    cells <- sample(1:N0, N)
  
    data.cs <- matrix(0, nrow=G0, ncol=N)
    
    for(k in 1:length(cells)){
      c.inx <- cells[k]
      
      # directly access to the sparse data...
      p.ks <- data@p[c.inx] # position start
      p.ke <- data@p[c.inx+1] # position end
      inx.k <- data@i[(p.ks+1):p.ke]+1 # indexes
      x.k <- data@x[(p.ks+1):p.ke] # data
      
      data.cs[inx.k, k] <- x.k # input the retrieved data
      
    }
    
    if(shuffle){
      rownames(data.cs) <- rownames(data)[sample(1:dim(data)[1], dim(data)[1])]
    }
    else{
      rownames(data.cs) <- rownames(data)
    }
    
    colnames(data.cs) <- paste0("cs_",m,"_",1:N)
    
    sobj.cs <- CreateSeuratObject(counts=data.cs)
    sobj.cs <- FindVariableFeatures(sobj.cs, selection.method = selection.method, nfeatures = nfeatures, verbose=verbose)
    # sobj.cs <- ScaleData(sobj.cs, verbose=verbose)
    
    vgenes <- VariableFeatures(object = sobj.cs)
    
    # GRN
    if(length(goi)==0){
      goi <- vgenes
    }
    vgenes.goi <- vgenes[vgenes%in%goi]
    
    subdata <- data.cs[vgenes.goi,]
    
    sd.mim <- build.mim(t(as.data.frame(subdata)), estimator = estimator, disc = disc, nbins=sqrt(ncol(subdata)))
    sd.aracne <- aracne(sd.mim, eps=10^-5)
    
    write.table(sd.aracne, paste0(directory,"cellstrap_", suffix, m,".txt"), quote=F, sep="\t")
    
  } 
  
}

# choice of estimator??
# "mi.empirical" : This estimator computes the entropy of the empirical probability distribution.
# "mi.mm" : This is the Miller-Madow asymptotic bias corrected empirical estimator.
# "mi.shrink" : This is a shrinkage estimate of the entropy of a Dirichlet probability distribution.
# "mi.sg" : This is the Schurmann-Grassberger estimate of the entropy of a Dirichlet probability
# distribution.
# "pearson" : This computes mutual information for normally distributed variable.
# "spearman" : This computes mutual information for normally distributed variable using Spearman’s
# correlation instead of Pearson’s correlation.
# "kendall" : This computes mutual information for normally distributed variable using Kendall’s
# correlation instead of Pearson’s correlation.

# usage
# annot <- read.table("c:/work/20230212_cellstrap/20230302_scrna_test_2/features_annot.tsv", row.names=1, header=F, fill=T, na.strings=NA, stringsAsFactors = F)
# goi <- read.table("c:/work/20230212_cellstrap/20230302_scrna_test_2/TF.geneids.txt", header=F, stringsAsFactors = F)
# 
# annot.goi <- annot[annot[,2]%in%unlist(goi),3]
# 
# cellstrap(data.kns, N=5000, M=30, goi=annot.goi, directory = "test/")
# 
# x <- read.table("cellstrap_1.txt", row.names=1, header=T, stringsAsFactors = F)

#===================================
# (3) aggregate the data; count TFs, count edges
# 11/14/2023: avoid loop...
# 11/14/2023: "nested" version; you can run the nested loop here
# R: number of nested iteration, 0 means no iteration

aggregateBCN <- function(directory, verbose=FALSE, R=100){
  files <- dir(directory)
  
  # first, check genes involved
  genes1 <- NULL; genes2 <- NULL
  for(f in files){
    x <- read.table(paste0(directory, f), row.names=1, header=T, stringsAsFactors = F)
    genes1 <- c(genes1, rownames(x))
    genes2 <- c(genes2, colnames(x))
  }
  genes.count <- table(genes1)
  genes1 <- sort(unique(genes1))
  genes2 <- sort(unique(genes2))
  # this step is quick
  
  # next, average matrices
  sd.aracne <- matrix(rep(0,length(genes1)*length(genes2)),ncol=length(genes1), nrow=length(genes2))
  count.aracne <- matrix(rep(0,length(genes1)*length(genes2)),ncol=length(genes1), nrow=length(genes2))
  rownames(sd.aracne) <- genes1; colnames(sd.aracne) <- genes2
  rownames(count.aracne) <- genes1; colnames(count.aracne) <- genes2
  
  min.aracne <- 1e6
  nf <- 0
  for(f in files){
    x <- read.table(paste0(directory, f), row.names=1, header=T, stringsAsFactors = F)
    min.aracne <- min(min(x[x>0]), min.aracne)
    genes.x1 <- rownames(x); genes.x2 <- colnames(x)
    
    sd.aracne[genes.x1, genes.x2] <- sd.aracne[genes.x1, genes.x2] + as.matrix(x)
    count.aracne[genes.x1, genes.x2] <- count.aracne[genes.x1, genes.x2] + (as.matrix(x)>0)
    
    # shuffle the genes if needed
    if(R>0){
      for(r in 1:R){
        shuffled.genes.x1 <- genes.x1[sample(1:length(genes.x1), length(genes.x1))]
        shuffled.genes.x2 <- genes.x2[sample(1:length(genes.x2), length(genes.x2))]        
        
        sd.aracne[shuffled.genes.x1, shuffled.genes.x2] <- sd.aracne[shuffled.genes.x1, shuffled.genes.x2] + as.matrix(x)
        count.aracne[shuffled.genes.x1, shuffled.genes.x2] <- count.aracne[shuffled.genes.x1, shuffled.genes.x2] + (as.matrix(x)>0)
      }
    }
    
    nf <- nf+1
    if(verbose){
      cat(nf)
    }
    
  }
  
  colnames(sd.aracne) <- rownames(sd.aracne)
  sd.aracne <- sd.aracne/nf
  if(R>0){
    sd.aracne <- sd.aracne/R
  }
  
  colnames(count.aracne) <- rownames(count.aracne)
  
  return(list(aracne=sd.aracne, node.count=genes.count, edge.count=count.aracne))
}

# usage
# agr.cells <- aggregate.cellstrap(directory="test/")
# 
# tiff("20230620_node-count_all.tiff")
# hist(agr.cells$node.count, breaks=30)
# dev.off()
# tiff("20230620_edge-count_all.tiff")
# hist(agr.cells$edge.count[agr.cells$edge.count>0], breaks=30)
# dev.off()


#===================================
# (4) p value?
# x: edge count
# alpha = X/MK*M0
# M0: number of resampling for x 
# X: shuffled edge count
# M: number of resampling for shuffled data X
# K: number of shuffling

computeEdgePvals <- function(x, alpha){
  
  genes.x1 <- rownames(x)
  genes.x2 <- colnames(x)
  genes.alpha <- rownames(alpha)
  
  pvals <- matrix(nrow=dim(x)[1], ncol=dim(x)[2])
  rownames(pvals) <- genes.x1
  colnames(pvals) <- genes.x2
  
  global.min.alpha <- min(alpha[alpha>0])
  
  for(gene1 in genes.x1){
    for(gene2 in genes.x2){
      
      # if nodes are not found in shuffled data, pval = minimum
      if(!(gene1 %in% genes.alpha)){
        pvals[gene1, gene2] <- 1-ppois(x[gene1, gene2], lambda=global.min.alpha)
        next
      }
      else if(!(gene2 %in% genes.alpha)){
        pvals[gene1, gene2] <- 1-ppois(x[gene1, gene2], lambda=global.min.alpha)
        next
      }
      
      # if the edge is not present in shuffled data, pval = minimum among edges connected to the gene
      if(alpha[gene1, gene2]==0){
        if(x[gene1, gene2]>0){
          min.alpha <- min((alpha[gene1,])[alpha[gene1,]>0])
          pvals[gene1, gene2] <- 1-ppois(x[gene1, gene2], lambda=min.alpha)
          next
        }
        # but, if the edge is not present in the current data, pval=1
        else{
          pvals[gene1, gene2] <- 1
          next
        }
      }
      
      pvals[gene1, gene2] <- 1-ppois(x[gene1, gene2], lambda=alpha[gene1, gene2])
      
    }
  }
  
  return(pvals)
  
}

#----------------------
# x: node count
# alpha = V/N*M
# M: number of resampling iteration
# V: number of variable expression genes selected 
# N: number of all the genes

computeNodePvals <- function(x, alpha){
  
  y <- x
  for(i in 1:length(x)){
    y[i] <- 1-ppois(x[i], lambda=alpha)
  }
  
  return(y)
  
}


#===================================
# (5) GRN graph
# cf C:/work/20230212_cellstrap/20231001_scrna_test_8_covid19/NatMed_2020_Liao/20231016_analysis_Liao-CoV_1.R
cleanGRN <- function(aracne.mtx, node.pvals, edge.pvals, node.cutoff=.0001, edge.cutoff=.0001){
  
  require(igraph)
  
  ar.graph0 <- graph.adjacency(aracne.mtx, mode="undirected", weighted=TRUE, diag=F)
  eweight0 <- log(igraph::E(ar.graph0)$weight)-floor(min(log(igraph::E(ar.graph0)$weight)))+1
  ar.graph0_fixedlayout <- layout_with_dh(ar.graph0)
  
  nodes.sig <- names(node.pvals)[node.pvals<node.cutoff] # take only significant nodes
  
  aracne.mtx.clean <- aracne.mtx
  aracne.mtx.clean[edge.pvals>edge.cutoff] <- 0 # remove non-significant edges
  aracne.mtx.clean <- aracne.mtx.clean[nodes.sig,nodes.sig] # remove non-significant nodes
  aracne.connectivity <- rowSums(aracne.mtx.clean>0)
  
  nodes.sig2 <- names(aracne.connectivity)[aracne.connectivity>0] # remove non-connected nodes
  aracne.mtx.clean <- aracne.mtx.clean[nodes.sig2,nodes.sig2]
  
  ar.graph <- graph.adjacency(aracne.mtx.clean, mode="undirected", weighted=TRUE, diag=F)
  ar.graph_fixedlayout <- layout_with_dh(ar.graph)
  
  eweight <- log(igraph::E(ar.graph)$weight)-floor(min(log(igraph::E(ar.graph)$weight)))+1
  
  ar.graph.list <- list(graph.all=list(graph=ar.graph0, weight=eweight0, layout=ar.graph0_fixedlayout),
                        graph.clean=list(graph=ar.graph, weight=eweight, layout=ar.graph_fixedlayout))
  
  return(ar.graph.list)
  
}



# usage
#
# pdf("20231016_cellstrap-graph0.pdf")
# plot(ar.graph$graph.clean$graph, layout=ar.graph$graph.clean$layout,
#      vertex.color="white", vertex.cex=0.5, vertex.size=10, vertex.label.cex=0.5,
#      edge.width=ar.graph$graph.clean$weight)
# dev.off()



#===================================
# writing MM to gz
writeMMgz <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    gzfile(file)
  )
  data.table::fwrite(
    x = summary(x),
    file = file,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}


