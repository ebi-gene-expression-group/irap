#!/usr/bin/env Rscript
#; -*- mode: R;-*-
# =========================================================
# Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================

###############################################################
suppressPackageStartupMessages(library("optparse"))

IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
if ( IRAP.DIR == "" ) {
  cat("ERROR: environment variable IRAP_DIR is not set\n")
  q(status=1)
}

# Load accessory functions

source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))
pdebug.enabled <- TRUE
#######################
usage <- "irap_sc3 --ifile file --out out_file_prefix [options]"
filenames <- c("ifile","lengths")

# specify our desired options in a list

option_list <- list(
    make_option(c("-i", "--ifile"), type="character", dest="ifile", help="TSV or MTX file name with the matrix with the counts per cell by gene/transcript."),
    make_option(c("-j", "--jfile"), type="character", dest="jfile", help="TSV or MTX file name with the matrix with the normalised counts per cell by gene/transcript."),
    make_option(c("-s", "--spikes"), type="character", dest="spikes", default=NULL, help="Optional GTF file with spike features"),
    make_option(c("--mtx"), action="store_true", default=FALSE,dest="is_mtx",help="The input file is in Matrix Market format. Default is TSV format."),
    make_option(c("--tsv"), action="store_false", default=FALSE,dest="is_mtx",help="The input file is in TSV format (default)."),
    make_option(c("-o", "--out"), type="character", default = NULL, help="Output file name prefix."),
    make_option(c("-d", "--dim"), type="numeric", default=3, help="Number of dimensions (default=3)."),
    make_option(c("-p", "--pca_components"), type="numeric", default = 30),
    make_option(c("-P", "--max-perplexity"), type="numeric",default=35,dest="max.perplexity",help="Maximum perplexity to consider."),
    make_option(c("--debug"), action="store_true", dest="debug", default=FALSE, help="Debug mode")
)

multiple.options <- NULL
mandatory <- c("ifile", "jfile", "out")

# Parse arguments

opt <- myParseArgs(
    usage = usage,
    option_list = option_list,
    filenames.exist = filenames,
    multiple.options = multiple.options,
    mandatory = mandatory
  )
pdebug.enabled <- opt$debug

pdebug.save.state("irap_scater_tsne", "p0")

pinfo("Loading ",opt$ifile)
if ( opt$is_mtx ) {
    raw <- mtx.load(opt$ifile)
    norm <- mtx.load(opt$jfile)
} else {
    raw <- as.matrix(quant.load(opt$ifile))
    norm <- as.matrix(quant.load(opt$jfile))
}
pinfo("Loading ",opt$ifile," done.")

# Load packages only after checking arguments and inputs

suppressPackageStartupMessages(library(scater))

sce <- SingleCellExperiment(
  assays = list(
    counts = raw,
    logcounts = norm
  )
)
rowData(sce)$feature_symbol <- rownames(sce)

# Identify the spikes (where present)

if ( ! is.null(opt$spikes) ){
  spikes <- load.gtf(opt$spikes)
  if (any(unique(spikes$seqid) %in% rownames(sce))){
    pinfo(paste("Marking", length(intersect(unique(spikes$seqid), rownames(sce))), "spikes"))
    isSpike(sce, "spikes") <- unique(spikes$seqid)
    
    # Not sure if Scater ignores the spikes for PCA etc, but let's take them out
    # anyway
    
    pinfo("Removing spikes")
    sce <- sce[! isSpike(sce), ]
    
  }else{
    pinfo(paste("Supplied spikes in", opt$spikes, "don't exist in expression matrix"))
    opt$spikes <- NULL
  }
}

# Run PCA 

pinfo("PCA...")

sce <- runPCA(sce, ncomponents = opt$pca_components)

# Run t-SNE for the range of perplexity values

perp.values <- c(1, seq(5,opt$max.perplexity,5))

for (p in perp.values) {
  pinfo("TSNE: perplexity ", p)
  sce <- runTSNE(sce, perplexity = p, rand_seed = 1000, exprs_values = 'logcounts', use_dimred = 'PCA', n_dimred = opt$pca_components, ncomponents = opt$dim )
  tsne_result <- data.frame(sce@reducedDims$TSNE)
  colnames(tsne_result) <- paste0('tSNE_', 1:ncol(tsne_result))
  tsne_result$Label <- rownames(tsne_result)
  
  ofile <- paste0(opt$out,"_tsne_perp_",p,".tsv")
  write.tsv(tsne_result, file=ofile)
}