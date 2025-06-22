#!/bin/sh

# compress and index vcf file
bgzip -c hgdp_tgp_sgdp_chr12_p.dated_287.vcf > file.vcf.gz
tabix -p vcf hgdp_tgp_sgdp_chr12_p.dated_287.vcf.gz

./bcftools/bcftools view -R SNPs_under_sel_3.txt -o ./human_data_SNPs_under_sel.vcf -O v ./hgdp_tgp_sgdp_chr12_p.dated_287.vcf.gz
