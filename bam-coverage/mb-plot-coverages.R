#==================#
#   HEADER START   #
#==================#
### Created: Jun 8, 2017
### Updates: Feb, 11 2019
### Author: Maciej Bak, Foivos Gypas
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Requirements: Gviz, rtracklayer, biomaRt, optparse, GenomicFeatures, GenomicRanges, GenomicAlignments
### R version used: 3.5.1
#==================#
### Description: Plots coverage of genes containing ROI
### Output: A directory with novel plots
#==================#
#    HEADER END    #
#==================#


# suppress warnings
options(warn=-1)

#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("Gviz"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("GenomicAlignments"))


#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
	make_option(c("--gtf"), action="store", type="character", default="", help="REQUIRED: GTF file with annotated transcripts", metavar="file"),
	make_option(c("--bed"), action="store", type="character", default="", help="REQUIRED: BED file with regions of interest", metavar="file"),
	make_option(c("--design_table"), action="store", type="character", default="", help="REQUIRED: TSV file with the information about alignment files", metavar="file"),
	make_option(c("--type"), action="store", type="character", default="", help="REQUIRED: Plot type: OVERLAY_ALL, OVERLAY_GROUPS, SEPARATE", metavar="string"),
	make_option(c("--output_dir"), action="store", type="character", default="", help="REQUIRED: Output directory", metavar="directory"),
	make_option(c("--help"), action="store_true", default=FALSE, help="Show this information and die"),
	make_option(c("--verbose"), action="store_true", default=FALSE, help="Be Verbose")
)
## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf [FILE] --bed [FILE] --design_table [FILE] --type [STRING] --output_dir [DIRECTORY]", option_list = option_list, add_help_option=FALSE, description="")
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if ( opt$gtf== "" || opt$bed=="" || opt$design_table=="" || opt$type=="" || opt$output_dir=="") {
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}
if (opt$type!="OVERLAY_ALL" && opt$type!="OVERLAY_GROUPS" && opt$type!="SEPARATE") {
	write("[ERROR] Invalid --type argument provided!\n\n", stderr())
	stop(print_help(opt_parser))
}
if ( opt$verbose ){
	options(warn=0)
}
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#

#---> START MESSAGE <---#
if ( opt$verbose  ) cat("Starting script", "'...\n\n", sep="")

#---> Configuration for non ucsc chromosomes (chomosome that start with chr) <---
options(ucscChromosomeNames=FALSE)

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose  ) cat("Reading annotation file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object
txdb <- makeTxDbFromGFF(opt$gtf)
geneTrack <- GeneRegionTrack(txdb,fontcolor="black",collapseTranscripts="meta",transcriptAnnotation="gene",col="#a9a9a9",fill="#a9a9a9")
names(geneTrack) = ""

#---> IMPORT bed regions <---#
# Print status message
if ( opt$verbose  ) cat("Reading bed file '", basename(opt$bed), "'...\n", sep="")
bed <- read.table(opt$bed,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

#---> Create output directory <---#
if ( opt$verbose    ) cat("Creating output directory'", basename(opt$output_dir), "'...\n", sep="")
dir.create(opt$output_dir, showWarnings = FALSE)

#---> Get the gene information <---#
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$gtf, format="gtf")
gr <- gr[values(gr)[["type"]] == "gene"]
gr <- gr[gr$gene_id %in% bed[["V9"]] ]

#---> Read the design table <---#
if ( opt$verbose  ) cat("Reading design table '", basename(opt$design_table), "'...\n", sep="")
design_table <- read.table(opt$design_table,header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="",comment.char="")

# Get the coverage counts
# The Alignmenttracks plot a smoothened version so this solution is not ideal:
# https://support.bioconductor.org/p/113298/#113363
coverages = apply(design_table, 1, function(x){
	return(coverage(readGAlignments(x["bam"]))[gr])
	})

if ( opt$verbose  ) cat("Generating plots")

# Separate tracks per every bam file
if (opt$type=="SEPARATE"){

	# plot regions one by one
	for (i in 1:nrow(bed)){

		svg(file.path(opt$output_dir,paste(i,"svg",sep=".")), width=10, height=4)

		#---> IMPORT alignment files <---#
		alignment_tracks = apply(design_table, 1, function(x){
			track = AlignmentsTrack(x["bam"], isPaired=FALSE, fill.coverage=x["color"], lwd.coverage="2", col.coverage=x["color"], alpha=1) 
			names(track) = x["sample"]
			return(track)
			})

		# Determine the max coverage in this region
		ylims = c(0,max(sapply(coverages, function(cvg){
			return(max(cvg[[toString(bed[i,"V1"])]]))
			})))

		# Put a frame around ROI
		ht = HighlightTrack(trackList = c(alignment_tracks,geneTrack), start=bed[i,"V7"], end=bed[i,"V8"], chromosome=bed[i,"V1"])
		displayPars(ht) <- list(col="black", lwd=1, fill="white")

		# Plot the tracks
		plotTracks(c(ht),
			chromosome=bed[i,"V1"],
			from=bed[i,"V2"],
			to=bed[i,"V3"],
			type="coverage",
			main=bed[i,"V4"],
			ylim=ylims,
			cex.main=1.0,
			sizes=c(rep(4,length(alignment_tracks)),1)
		)

		dev.off()
	}

# Overlay all tracks as coverage contours
} else if (opt$type=="OVERLAY_ALL"){

	# plot regions one by one
	for (i in 1:nrow(bed)){

		svg(file.path(opt$output_dir,paste(i,"svg",sep=".")), width=10, height=4)

		#---> IMPORT alignment files <---#
		alignment_tracks = apply(design_table, 1, function(x){
			track = AlignmentsTrack(x["bam"], isPaired=FALSE, fill.coverage="transparent", lwd.coverage="2", col.coverage=x["color"], alpha=1) 
			names(track) = x["sample"]
			return(track)
			})

		# Determine the max coverage in this region
		ylims = c(0,max(sapply(coverages, function(cvg){
			return(max(cvg[[toString(bed[i,"V1"])]]))
			})))

		# overlay all alignment tracks
		ot <- OverlayTrack(trackList = alignment_tracks)

		# Put a frame around ROI
		ht = HighlightTrack(trackList = c(ot,geneTrack), start=bed[i,"V7"], end=bed[i,"V8"], chromosome=bed[i,"V1"])
		displayPars(ht) <- list(col="black", lwd=1, fill="white")

		# Plot the tracks
		plotTracks(c(ht),
			chromosome=bed[i,"V1"],
			from=bed[i,"V2"],
			to=bed[i,"V3"],
			type="coverage",
			main=bed[i,"V4"],
			ylim=ylims,
			cex.main=1.0,
			sizes=c(4*length(alignment_tracks),1)
		)

		dev.off()
	}

# Overlay coverages from the same group, plot groups as separate tracks
} else {
	stopifnot(opt$type=="OVERLAY_GROUPS")

	# plot regions one by one
	for (i in 1:nrow(bed)){

		svg(file.path(opt$output_dir,paste(i,"svg",sep=".")), width=10, height=4)

		# Determine the max coverage in this region
		ylims = c(0,max(sapply(coverages, function(cvg){
			return(max(cvg[[toString(bed[i,"V1"])]]))
			})))

		# create separate overlaid tracks per every group
		#---> IMPORT alignment files <---#
		split_df = split.data.frame(design_table,design_table$group)
		group_list = lapply(split_df,"[[","bam")
		ots = lapply(names(group_list), function(x){
			alignment_tracks = lapply(group_list[[x]], function(p){
				color = c(design_table[grepl(p,design_table$bam),]["color"])$color
				AlignmentsTrack(p, isPaired=FALSE, fill.coverage="transparent", lwd.coverage="2", col.coverage=color, alpha=1)
				})
			track = OverlayTrack(trackList = alignment_tracks)
			names(track) = x
			return(track)
			})

		# Put a frame around ROI
		ht = HighlightTrack(trackList = c(ots,geneTrack), start=bed[i,"V7"], end=bed[i,"V8"], chromosome=bed[i,"V1"])
		displayPars(ht) <- list(col="black", lwd=1, fill="white")

		# Plot the tracks
		plotTracks(c(ht),
			chromosome=bed[i,"V1"],
			from=bed[i,"V2"],
			to=bed[i,"V3"],
			type="coverage",
			main=bed[i,"V4"],
			ylim=ylims,
			cex.main=1.0,
			sizes=c(rep(16/length(group_list),length(group_list)),1)
		)

		dev.off()
	}
}

#================#
#   MAIN END     #
#================#
