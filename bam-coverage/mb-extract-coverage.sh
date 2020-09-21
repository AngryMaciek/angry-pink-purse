#
# Maciek Bak
# 14-APR-2019
#
# [USAGE]
# bash mb-extract-coverage.sh bam=<BAM> bed=<BED> genome=<FASTA> outdir=<PATH>
#
# [REQUIREMENTS]
# samtools=1.9
# mosdepth=0.2.4
# bedtools=2.26.0
#
# [DESCRIPTION]
# Extract genomic coverages given in bed file from bam file.
# For the interpretation of the output files look here:
# https://github.com/brentp/mosdepth
# =============================================================================

# evaluate the command-line arguments (based on order)
eval $1
eval $2
eval $3
eval $4

# check if bed and bam file exist
if [ ! -e $bam ]; then
    echo "### BAM file not found!"
    exit 1
fi

if [ ! -e $bed ]; then
    echo "### BED file not found!"
    exit 1
fi

# create the output directory, if exists - exit
mkdir $outdir
if [ $? -ne 0 ] ; then
    echo "### Output directory exists!"
    exit 1
fi

# sort the bam file by position
samtools sort $bam 1> $outdir/sample_sorted.bam 2> /dev/null;

# create the index for sorted bam
samtools index $outdir/sample_sorted.bam 1> $outdir/sample_sorted.bam.bai;

# sort the bed file by coordinates
sort -k1,1 -k2,2n $bed 1> $outdir/regions.bed;

# calculate the genomic coverage
mosdepth --by $outdir/regions.bed $outdir/sample $outdir/sample_sorted.bam;

# generate fasta index for the genome
# reformat the index to get the bedtools format genome file
samtools faidx $genome;
mv $genome.fai $outdir/genome.fai;
awk -v OFS='\t' {'print $1,$2'} $outdir/genome.fai > $outdir/genome_file.txt;
rm -f $outdir/genome.fai;

# intersect bed and bam files -> extract ROI from bam coverages
bedtools intersect -sorted \
-b $outdir/regions.bed -a $outdir/sample.per-base.bed.gz \
-g $outdir/genome_file.txt -wb > $outdir/intersect_coverage.bed;

# reformat the output to a regular bed
cut -f1,2,3,4,8,10 cov/intersect_coverage.bed | \
awk -v OFS='\t' ' { t = $4; $4 = $5; $5 = t; print; } ' \
> $outdir/ROI_coverage.bed
