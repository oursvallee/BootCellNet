# What's here

This is the repository for BootCellNet, a procedure for scRNA-seq analysis, containing several functions written in R.

bootcellnet2.R: contains all the functions

tutorial_files.tar.gz: containing files necessary for test BootCellNet

# Inquiries

For more details, please see

BootCellNet, a resampling-based procedure, promotes unsupervised identification of cell populations via robust inference of gene regulatory networks.

https://www.biorxiv.org/content/10.1101/2024.02.06.579236v2

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

Here's a tutorial for those who want to test BootCellNet.

## Install

Download tutorial_files.tar.gz, and expand all the files into one directory. Put bootcellnet2.R into the same directory.

PBMC_test.R: A kind of tutorial. You can run BootCellNet on PBMC data provided by 10x genomics to see how it works.

Also, please download the required data from 10x.

wget https://cf.10xgenomics.com/samples/cell-exp/3.1.0/manual_5k_pbmc_NGSC3_ch5/manual_5k_pbmc_NGSC3_ch5_filtered_feature_bc_matrix.tar.gz
gunzip manual_5k_pbmc_NGSC3_ch5_filtered_feature_bc_matrix.tar.gz

### Depending on the version of Seurat you are using, the UMAP image of the data may be different from what is shown in the BootCellNet paper.###
For convenience, I put RDS files and .RData files for both versions of Seurat. All are included in tutorial_files.tar.gz

Please check the comments in PBMC_test.R for detail.





