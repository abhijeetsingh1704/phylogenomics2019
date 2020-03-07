#!/bin/bash

###
# This script will summarize output files for gene-wise
# and site-wise log likelihood values and calculate
# gene-wise and site-wise log likelihood scores
# according to Shen et al. 2017
#
# A pdf of this manuscript is available here
# https://s3.amazonaws.com/vu-wp0/wp-content/uploads/sites/191/2017/04/10171603/2017_Shen_etal_NEE.pdf 
#
# The output from the script is ordered in the following:
# col 1: gene partition
# col 2: gene log likelihood for the reference topology
# col 3: gene log likelihood for the alternative topology
# col 4: the difference in log likelihood values
# col 5: the tree that is supported, ref or alt
#
# usage: bash create_GLS_summary.sh arg1 arg2 arg3 arg4 arg5
# arg 1: partition file
# arg 2: gene-wise log likelihood ref
# arg 3: gene-wise log likelihood alt
# arg 4: site-wise log likelihood ref
# arg 5: site-wise log likelihood alt
#
###

# save start time
start=$(date +%s.%N)

# help message
if [ "$1" == "-h" ]; then
  echo -e "\nThis script will summarize output files for gene-wise"
  echo "and site-wise log likelihood values and calculate"
  echo "gene-wise and site-wise log likelihood scores"
  echo -e "according to Shen et al. 2017\n"
  echo "The output from the script is ordered in the following:"
  echo "col 1: gene partition"
  echo "col 2: gene log likelihood for the reference topology"
  echo "col 3: gene log likelihood for the alternative topology"
  echo "col 4: the difference in log likelihood values"
  echo -e "col 5: the tree that is supported, ref or alt\n"
  echo "usage: bash create_GLS_summary.sh arg1 arg2 arg3 arg4 arg5"
  echo "arg 1: partition file"
  echo "arg 2: gene-wise log likelihood ref"
  echo "arg 3: gene-wise log likelihood alt"
  echo "arg 4: site-wise log likelihood ref"
  echo -e "arg 5: site-wise log likelihood alt\n"
  exit 0
fi

# read in arguments
echo "Reading arguments"
echo -e "..."
PART=$1
GENEref=$2
GENEalt=$3
SITEref=$4
SITEalt=$5
echo "Partition file: $PART"
echo "Gene-wise log likelihood ref: $GENEref"
echo "Gene-wise log likelihood alt: $GENEalt"
echo "Site-wise log likelihood ref: $SITEref"
echo "Site-wise log likelihood alt: $SITEalt"
echo -e "Complete!\n\n"

# create gene-wise log likelihood summary file
NUMGENES=$(cat $GENEref | tr " " "\n" |  grep "-" | wc -l)
echo "Creating the gene-wise log likelihood summary file for $NUMGENES genes"
echo -e "..."
paste <(cat $PART | awk '{print $2}') \
  <(cat $GENEref | tr " " "\n" |  grep "-") \
  <(cat $GENEalt | tr " " "\n" |  grep "-") | \
  awk -v OFS='\t' '{print $0, $2-$3}' | \
  awk -v OFS='\t' '{if ($NF<0) print $0"\talt"; else if ($NF==0) print $0"\tneither"; else print $0"\tref"}' \
  > gene-wise_logli_scores.txt
echo -e "Complete!\n\n"

# create site-wise log likelihood summary file
NUMSITES=$(cat $SITEref | tr " " "\n" | grep "-" | wc -l)
echo "Creating the site-wise log likelihood summary file for $NUMSITES sites"
echo -e  "..."
paste <(seq 1 $NUMSITES) \
  <(cat $SITEref | tr " " "\n" |  grep "-") \
  <(cat $SITEalt | tr " " "\n" |  grep "-") | \
  awk -v OFS='\t' '{print $0, $2-$3}' | \
  awk -v OFS='\t' '{if ($NF<0) print $0"\talt"; else if ($NF==0) print $0"\tneither"; else print $0"\tref"}' \
  > site-wise_logli_scores.txt
echo -e "Complete!\n\n"

# report execution time 
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time: %.6f seconds\n" $dur
