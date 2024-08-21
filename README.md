# What's here

This is the repository for BootCellNet, a procedure for scRNA-seq analysis, containing several functions written in R.

bootcellnet2.R: contains all the functions

tutorial_files.tar.gz: containing files necessary for test BootCellNet

# Inquiries

For more details, please see

BootCellNet, a resampling-based procedure, promotes unsupervised identification of cell populations via robust inference of gene regulatory networks.

https://www.biorxiv.org/content/10.1101/2024.02.06.579236v1

or contact the author.

You can ask any question you may have from the anonymous Google Forms.

https://docs.google.com/forms/d/e/1FAIpQLScGG2CjumDr9DiXzl7ZBrMqxk2BuTc4fObHIjmf3cxORgwYdA/viewform?usp=sf_link

Feel free to ask anything. The answer will be posted here.

# How to use

## Prerequisite

It can be run on R. Additionally you need the following packages:

Seurat, Matrix, igraph, scales, ggplot2, gplots, minet, pvclust, lpSolve

## Install

Just download the "bootcellnet2.R", read it by source("YOUR_PATH/bootcellnet.R"), and perform any analysis you want to do.

# Tutorial

A tutorial has been provided for those who want to test BootCellNet.

## Install

Download tutorial_files.tar.gz, and expand all the files into one directory. Put bootcellnet2.R into the same directory.

Also, please download the required data from 10x before start:

wget https://cf.10xgenomics.com/samples/cell-exp/3.1.0/manual_5k_pbmc_NGSC3_ch5/manual_5k_pbmc_NGSC3_ch5_filtered_feature_bc_matrix.tar.gz
gunzip manual_5k_pbmc_NGSC3_ch5_filtered_feature_bc_matrix.tar.gz

### What's inside

PBMC_tutorial.R: A kind of tutorial. You can run BootCellNet on PBMC data provided by 10x genomics to see how it works.

MORF_genelist_geneids.txt: A list of GeneIDs for human transcription factors, taken from Jeong et al, Cell, 2024

pbmc_annot_geneid.txt: A list of annotation data for the 10x data.

#### Note
Depending on the version of Seurat you are using, the UMAP image of the data may be different from what is shown in the BootCellNet paper. Thus I included these files for the calculated results.

pbmc_kns_[old|new].mtx.gz: kNN-smoothed data, obtained from data processed by Seurat ver 4 or 5, respectively.

PBMC_tutorial_[old|new].RData: Results of calculation obtained by Seurat ver 4 or 5, respectively. This does not include 10x data.

pbmc_bcn_[old|new]_[|shuffle]: Folders containing results of bootstrap calculation of GRN, for the reconstruction (100 files) or the shuffling (1000 files), obtained from data processed by Seurat ver 4 or 5, respectively.

Please check the comments in PBMC_test.R for detail.

These files are not provided here, because of the space limit.

PBMC_tutorial_[old|new].rds: RDS file, containing 10x data, PCA and UMAP coordinates, and the results of clustering, calculated by Seurat ver 4 or 5, respectively.

Therefore, I made them available from the link below:

old: https://drive.google.com/file/d/1IuueU-eRHicJhNC5CoupmO4NiC6hQpnN/view?usp=sharing

new: https://drive.google.com/file/d/1MsHMxvSxyCU-DWbbJ5LGZCNsNLEXKltJ/view?usp=sharing




