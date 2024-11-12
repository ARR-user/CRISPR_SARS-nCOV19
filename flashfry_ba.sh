#!/bin/bash

#---------------------------- Arrays for chromosomes and genome IDs---------------------------------------------
chromosomes=("chr1" "chr2" )  
genomes=("HG00142" "HG00178" "HG00237" "HG01521" "HG03129" "NA19175" "NA19462")  

# Create output folder 
mkdir -p output

# Loop through chromosomes and genome IDs
for chr in "${chromosomes[@]}"; do
  for genome in "${genomes[@]}"; do
    # Set database and reference file paths
    database="${chr}_${genome}_cas9ngg_database"
    reference="${chr}_genome/${genome}_${chr}.fasta"

    # Step 1: Index
    echo "Indexing for $chr and $genome..."
    java -Xmx4g -jar FlashFry-assembly-1.15.jar \
      index \
      --tmpLocation ./tmp \
      --database "$database" \
      --reference "$reference" \
      --enzyme spcas9ngg
    echo "fin-1"  

    # Step 2: Discover
    echo "Discovering targets for $chr and $genome..."
    java -Xmx4g -jar FlashFry-assembly-1.15.jar \
      discover \
      --database "$database" \
      --fasta LPL_first_exon_HPSE_NC_000004.12.fasta \
      --output "output/${genome}_${chr}_LPL.output"
     echo "fin-2" 

    # Step 3: Score
    echo "Scoring targets for $chr and $genome..."
    java -Xmx4g -jar FlashFry-assembly-1.15.jar \
      score \
      --input "${genome}_${chr}_LPL.output" \
      --output "output/${genome}_${chr}_LPL.output.scored" \
      --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
      --database "$database"
       echo "fin-3" 

  done
done

echo "All tasks completed."
#  output text file in respective folders 
