#!/usr/bin/env Rscript
libloc = file.path("~", "R","library", fsep = .Platform$file.sep)
req = c("igraph", "slam", "BiRewire", "optparse", "tsne")
miss <- setdiff(req,  installed.packages(lib.loc = libloc)[,"Package"])
if (length(miss)) {
	source("http://bioconductor.org/biocLite.R")
	biocLite(lib=libloc, lib.loc=libloc)
	biocLite("BiRewire", lib=libloc, lib.loc=libloc, dependencies=TRUE)
	install.packages(c("igraph","slam", "optparse", "tsne"), lib=libloc) #removed lib.loc = libloc
}
library(BiRewire, lib.loc = libloc)
library(igraph)
library(slam)
library(optparse, lib.loc = libloc)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="path to incidence matrix file", metavar="character")  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

with (opt,{
	print (file)
	data <- read.table(file, sep="", colClasses='character')
	dat <-sapply(data, function(e) {
	  splitter <- strsplit(e,"")
	})
	dat <- sapply(dat, function (e) {
	as.numeric (unlist(e))
	})
	dat <- t(dat)
	sink(file=file)
	m2<-birewire.rewire.bipartite(dat,verbose=FALSE)
	subs_on_node <-apply(m2, 1, function(e){
	  which (e>0)
	})
	sapply(subs_on_node, function(e){
	  cat (paste(e, collapse=","), "\n") 
	})
	sink()
})