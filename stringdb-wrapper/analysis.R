#==================#
#   HEADER START   #
#==================#
### Created: Nov 16, 2018
### Author: Maciej Bak
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Description: StringDB analysis of a set of genes
#==================#
#    HEADER END    #
#==================#

# Based on:
# https://bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf

########################
### LOAD PACKAGES    ###
########################
library("optparse")
library("STRINGdb")
library("biomaRt")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
	make_option(c("--output_dir"), action="store", type="character", default="", help="REQUIRED: Output directory", metavar="directory"),
	make_option(c("--genes"), action="store", type="character", default="", help="REQUIRED: List of gene IDs"),
	make_option(c("--species"), action="store", type="character", default="", help="REQUIRED: species ID (hsa or mmu)"),
	make_option(c("--help"), action="store_true", default=FALSE, help="Show this information and die"),
	make_option(c("--verbose"), action="store_true", default=FALSE, help="Be Verbose")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --genes [FILE] --output_dir [DIRECTORY]", option_list = option_list, add_help_option=FALSE, description="")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if (opt$genes=="" || opt$output_dir=="") {
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}
if ( opt$verbose ){
	options(warn=0)
}

# create the output directory
dir.create(opt$output_dir)

# select the database for human/mouse
if (opt$species == "hsa"){

	# define biomart object : ENSEMBL names ENCODING
	mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

	# initiate the reference class object
	string_db = STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory=opt$output_dir)
	# string_db <- STRINGdb$new( score_threshold=0, backgroundV = backgroundV ) #PROVIDE BACKGROUND

} else if (opt$species == "mmu"){

	# define biomart object : ENSEMBL names ENCODING
	mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

	# initiate the reference class object
	string_db = STRINGdb$new(version="10", species=10090, score_threshold=0, input_directory=opt$output_dir)
	# string_db <- STRINGdb$new( score_threshold=0, backgroundV = backgroundV ) #PROVIDE BACKGROUND

}

# get the STRING IDs for our genes
genes = read.table(opt$genes)
names(genes) = c("my_IDs")
genes = string_db$map(genes,"my_IDs",takeFirst=TRUE,removeUnmappedRows=FALSE,quiet=FALSE)

# filter out genes with no products found in db
genes = genes[complete.cases(genes), ]

score_treshold=0

# generate the network plots

pdf(file=paste(opt$output_dir,"network.pdf",sep="/"))
par(mar=c(0,0,0,0))
string_db$plot_network(genes$STRING_id,required_score=score_treshold,add_link=FALSE,add_summary=FALSE)
dev.off()

pdf(file=paste(opt$output_dir,"network_annotated.pdf",sep="/"))
par(mar=c(4,4,4,4))
string_db$plot_network(genes$STRING_id,required_score=score_treshold,add_link=TRUE,add_summary=TRUE)
dev.off()

png(file=paste(opt$output_dir,"network.png",sep="/"), width=1024, height=1024)
par(mar=c(0,0,0,0))
string_db$plot_network(genes$STRING_id,required_score=score_treshold,add_link=FALSE,add_summary=FALSE)
dev.off()

png(file=paste(opt$output_dir,"network_annotated.png",sep="/"), width=1024, height=1024)
par(mar=c(4,4,4,4))
string_db$plot_network(genes$STRING_id,required_score=score_treshold,add_link=TRUE,add_summary=TRUE)
dev.off()

enrichmentGO_Process <- string_db$get_enrichment(genes$STRING_id, category = "Process", methodMT = "fdr", iea = TRUE)
enrichmentGO_Component <- string_db$get_enrichment(genes$STRING_id, category = "Component", methodMT = "fdr", iea = TRUE)
enrichmentGO_Function <- string_db$get_enrichment(genes$STRING_id, category = "Function", methodMT = "fdr", iea = TRUE)
enrichmentKEGG <- string_db$get_enrichment(genes$STRING_id, category = "KEGG", methodMT = "fdr", iea = TRUE)

write.table(enrichmentGO_Process, file=paste(opt$output_dir,"GO_Processess.tsv",sep="/"), sep="\t", quote=FALSE)
write.table(enrichmentGO_Component, file=paste(opt$output_dir,"GO_Component.tsv",sep="/"), sep="\t", quote=FALSE)
write.table(enrichmentGO_Function, file=paste(opt$output_dir,"GO_Function.tsv",sep="/"), sep="\t", quote=FALSE)
write.table(enrichmentKEGG, file=paste(opt$output_dir,"KEGG.tsv",sep="/"), sep="\t", quote=FALSE)

# Get the clusters of interacting proteins:

clusters_dir = paste(opt$output_dir,"clusters",sep="/")
dir.create(clusters_dir)

clustersList = string_db$get_clusters(genes$STRING_id)
# plot first n_clusters
n_clusters = 5

for(i in seq(1:n_clusters)){

  # save the genes of the cluster
  cluster_proteins = sapply(clustersList[[i]], function(x) strsplit(x,".",fixed=TRUE)[[1]][2])
  gene_id = getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                  filters = "ensembl_peptide_id", values = cluster_proteins, mart = mart)
  gene_id = gene_id$ensembl_gene_id
  fname = paste(paste("cluster",as.character(i),sep="_"),"txt",sep=".")
  write(gene_id, file = paste(clusters_dir,fname,sep="/"), ncolumns = 1, append = FALSE)
  
  # get all the interactors of the proteins of the cluster...
  neighbours = string_db$get_neighbors( clustersList[[i]] )
  neighbours = sapply(neighbours, function(x) strsplit(x,".",fixed=TRUE)[[1]][2])
  # ...and save all their gene IDs
  gene_id = getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                 	filters = "ensembl_peptide_id", values = neighbours, mart = mart)
  gene_id = gene_id$ensembl_gene_id
  fname = paste(paste("neighbours_for_cluster",as.character(i),sep="_"),"txt",sep=".")
  write(gene_id, file = paste(clusters_dir,fname,sep="/"), ncolumns = 1, append = FALSE)

	# get detailed info about the interactions
  interactions = string_db$get_interactions(clustersList[[i]])
  interactions$to = sapply(interactions$to, function(x) strsplit(x,".",fixed=TRUE)[[1]][2])
  protein_gene_dict = getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                            filters = "ensembl_peptide_id", values = interactions$to, mart = mart)
  colnames(protein_gene_dict) = c("to_geneID","ensembl_peptide_id")
  interactions = merge(x=interactions, y=protein_gene_dict, by.x="to", by.y="ensembl_peptide_id")  
  interactions$from = sapply(interactions$from, function(x) strsplit(x,".",fixed=TRUE)[[1]][2])
  protein_gene_dict = getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                 	filters = "ensembl_peptide_id", values = interactions$from, mart = mart)
  colnames(protein_gene_dict) = c("from_geneID","ensembl_peptide_id")
  interactions = merge(x=interactions, y=protein_gene_dict, by.x="from", by.y="ensembl_peptide_id")

  fname = paste(paste("interactions_for_cluster",as.character(i),sep="_"),"txt",sep=".")
  write.table(interactions, file=paste(clusters_dir,fname,sep="/"), append=FALSE, quote=FALSE, sep="\t")

  # make network plots for a cluster
  
  fname = paste(paste("cluster",as.character(i),sep="_"),"pdf",sep=".")
  pdf(file=paste(clusters_dir,fname,sep="/"))
  par(mar=c(0,0,0,0))
  string_db$plot_network(clustersList[[i]],required_score=score_treshold,add_link=FALSE,add_summary=FALSE)
  dev.off()	
	
  fname = paste(paste("cluster",as.character(i),sep="_"),"png",sep=".")
  png(file=paste(clusters_dir,fname,sep="/"), width=1024, height=1024)
  par(mar=c(0,0,0,0))
  string_db$plot_network(clustersList[[i]],required_score=score_treshold,add_link=FALSE,add_summary=FALSE)
  dev.off()	
	
}
