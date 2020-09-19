###############################################################################
#
#   A small bash script to extract only canonical chromosome entries (+MT)
#   from ENSEMBL genomic annotation in GTF format.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 26-11-2019
#   LICENSE: GPL_v3.0
#   USAGE:
#   bash filter.sh -s {hsa|mmu} -gtf {} -o {}
#
###############################################################################

###############################################################################
# parse command line arguments
###############################################################################

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -s|--species)
        SPECIES="$2"
        shift # past argument
        shift # past value
        ;;
    -gtf|--genomic_annotation)
        GTF="$2"
        shift # past argument
        shift # past value
        ;;
    -o|--output_file)
        OUTFILE="$2"
        shift # past argument
        shift # past value
        ;;
    *) # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

###############################################################################
# MAIN
###############################################################################

# check species option value
if [ ! ${SPECIES} == "hsa" ] && [ ! ${SPECIES} == "mmu" ]; then
    echo $"Invalid species option. Currently supported are: 'hsa', 'mmu'"
    exit 1
fi

# write GTF headers
awk '$1 ~ /^#/ {print $0;next}' ${GTF} > ${OUTFILE}

# add canonical chromosomes
case $SPECIES in
    hsa) # Homo sapiens
        awk '$1 {if ($1 == "1") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "2") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "3") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "4") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "5") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "6") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "7") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "8") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "9") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "10") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "11") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "12") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "13") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "14") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "15") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "16") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "17") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "18") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "19") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "20") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "21") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "22") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "X") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "Y") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "MT") print}' ${GTF} >> ${OUTFILE}
        ;;
    mmu) # Mus musculus
        awk '$1 {if ($1 == "1") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "2") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "3") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "4") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "5") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "6") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "7") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "8") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "9") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "10") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "11") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "12") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "13") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "14") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "15") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "16") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "17") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "18") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "19") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "X") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "Y") print}' ${GTF} >> ${OUTFILE}
        awk '$1 {if ($1 == "MT") print}' ${GTF} >> ${OUTFILE}
        ;;
esac
