###############################################################################
#
#   Bash script to convert ENSEMBL genomic annotation from GTF format to BED.
#
#   This script requires package BEDOPS loaded in the current environment.
#   https://bedops.readthedocs.io/en/latest/index.html
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 12-02-2020
#   LICENSE: GPL_v3.0
#   USAGE: mb_convert_ENSEMBL_GTF_2_BED.sh -gtf {GTF} -bed {BED}
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
    -gtf)
        GTF="$2"
        shift # past argument
        shift # past value
        ;;
    -bed)
        BED="$2"
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

grep -v "^#\!" ${GTF} \
| awk \
'{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' \
| gtf2bed - \
1> ${BED}
