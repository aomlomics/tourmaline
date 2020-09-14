library(msa)
library(odseq)
myAlignment = readDNAMultipleAlignment(
	filepath = snakemake@input[[1]],
	format = "fasta")
myOutliers = odseq(
	myAlignment,
	distance_metric = snakemake@params[["metric"]],
	B = snakemake@params[["replicates"]],
	threshold = snakemake@params[["threshold"]])
write.table(
	myOutliers,
	snakemake@output[[1]],
	sep = ",",
	col.names = FALSE,
	quote = FALSE)
