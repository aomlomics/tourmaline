#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(msa)
library(odseq)
myAlignment = readDNAMultipleAlignment(
	filepath = args[1],
	format = "fasta")
myOutliers = odseq(
	myAlignment,
	distance_metric = args[2],
	B = as.numeric(args[3]),
	threshold = as.numeric(args[4]))
write.table(
	myOutliers,
	args[5],
	sep = "\t",
	col.names = FALSE,
	quote = FALSE)
