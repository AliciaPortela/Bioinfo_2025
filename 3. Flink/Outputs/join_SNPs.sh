#!/bin/sh

# join all SNPs under balancing selection to filter original VCF with bcftools
awk '{print $2 "\t" $4}' coefSelbalancing* | sort -k2,2n | uniq -f1 > all_balanced_loci.txt

# join all SNPs under divergent selection to filter original VCF with bcftools
awk '{print $2 "\t" $4}' coefSeldivergent* | sort -k2,2n | uniq -f1 > all_divergent_loci.txt
