#!/bin/bash

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 6 ]; then
cat << EOF
Usage: sh ${SCRIPT} vcf_dir pop_file pop_prefix out_dir cores downsample_num.  Parameters must be passed in this order

All arguments are required and must follow correct order.

vcf_dir: input directory containing VCFs to be subset by population. 

pop_file: file that links sample IDs in the VCF to populations 

pop_prefix: population prefix in the pop_file to be subset

out_dir: Directory where the subsetted VCFs will be output

cores: Number of processors to use for parallelization

downsample_num: Number of individuals to downsample each population to.  Set to -9 to shut off downsampling

EOF
  exit 1
fi

# check that input directory exists and is a directory
if ! [ -d "$1" ]; then
  echo "$1 is not a directory" >&2
  exit 1
fi
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi

# Check that population file exists
if ! [ -e "$2" ]; then
  echo "$2 not found" >&2
  exit 1
fi

# make output directory if it doesn't already exist
if ! [ -e "$4" ]; then
  mkdir -p $4
fi

module load parallel

VCF_DIR="$1"
PFILE="$2"
POP="$3"
OUT="$4"
THREADS="$5"
DS="$6"

#New file containing just the sample IDs for the population of interest
SFILE=$(mktemp ${OUT}/SAMPS.XXXXXXXXX)

#Retrieve sample IDs for this population
if [ "$DS" -ne "-9" ]; then 
  grep -w ${POP} ${PFILE} | awk '{print $1}' | shuf -n $DS > ${SFILE}
else
  grep -w ${POP} ${PFILE} | awk '{print $1}' > ${SFILE}
fi


POP_SUBSET() {
  module load bcftools/1.9
  module load htslib/1.9
  source activate py27

  IN=$1
  OUT=$2
  SFILE=$3
  POP=$4

  BS=$(basename ${IN})

  DUPS=$(mktemp ${OUT}/dups.XXXXXXXXX)

  zgrep -v "#" $IN | awk '{print $3}' | sort | uniq -d > $DUPS
  bcftools view --force-samples -S ${SFILE} -Oz -o ${OUT}/${BS%.vcf.gz}.${POP}.vcf.gz ${IN} &&
  tabix -p vcf ${OUT}/${BS%.vcf.gz}.${POP}.vcf.gz;
  #Make additional file that is filtered to only contain SNPs where ancestral state has been labelled
  zcat ${OUT}/${BS%.vcf.gz}.${POP}.vcf.gz | vawk --header '{if(I$AA=="a|||" || I$AA=="A|||" || I$AA=="g|||" || I$AA=="G|||" || I$AA=="t|||" || I$AA=="T|||" || I$AA=="c|||" || I$AA=="C|||") print $0}' | grep -v -f $DUPS | bgzip > ${OUT}/${BS%.vcf.gz}.${POP}.AAflt.vcf.gz;
  mv ${OUT}/${BS%.vcf.gz}.${POP}.AAflt.vcf.gz ${OUT}/${BS%.vcf.gz}.${POP}.vcf.gz
  tabix -p vcf ${OUT}/${BS%.vcf.gz}.${POP}.vcf.gz
  rm $DUPS

}

export -f POP_SUBSET

parallel -j $THREADS POP_SUBSET {} $OUT $SFILE $POP ::: ${VCF_DIR}/*vcf.gz

rm $SFILE






