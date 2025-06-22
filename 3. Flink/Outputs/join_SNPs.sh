#!/bin/sh

# join all SNPs under balancing selection to filter original VCF with bcftools
awk '{print $2 "\t" $4}' coefSelbalancing* > all_balanced_loci.txt

# join all SNPs under divergent selection to filter original VCF with bcftools
awk '{print $2 "\t" $4}' coefSeldivergent* > all_divergent_loci.txt
