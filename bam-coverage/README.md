# Toolset to extract and plot genomic coverages from bam files
*Maciej Bak  
maciej.bak@unibas.ch  
7.02.2019*

Extracting genomic coverage is not a straightforward task:

https://www.biostars.org/p/60841/  
https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html  
http://bioconductor.org/packages/release/bioc/html/bamsignals.html  
https://biostar.usegalaxy.org/p/18245/

## mb-plot-coverages.R

R script that takes bam files and plots the coverage profiles.

(Decription inside the script)

Input:

* design table with path to the bam files, sample names and colors:

  ```
  sample	bam	group	color
  [sample_name]	[path_to_bam_file]	[group_name]	[hex_color]
  ```

* extended bed file with the regions to plot and region to highlight

  ```
  [chr]	[plot_region_start]	[plot_region_end]	[region_name]	[score(unused)]	[strand]	[highlight_region_start]	[highlight_region_end]	[ENSEMBL_gene_ID]
  ```

* annotation file in the gtf format

Sample call:

```
Rscript mb-plot-coverages.R \
--gtf ~/RESOURCES/Homo_sapiens.GRCh38.87.gtf \
--bed regions.bed \
--design_table design_table.tsv \
--type SEPARATE \
--output_dir .
```

## mb-extract-coverage.sh

Bash script that takes a sorted bam file and a standard bed file and produces the position-wise coverage in the regions of interest.

(Decription inside the script)

Sample call:

```bash
bash mb-extract-coverage.sh \
bam=star/WT1/STAR_Aligned.out.sorted.bam \
bed=regions.bed \
outdir=cov \
genome=Homo_sapiens.GRCh38.dna.primary_assembly.fa
```
