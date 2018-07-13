#!/usr/bin/evn
ad='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/test_ad'

bcftools view -h batch1_75_AD_SNPs_GQ50_AD08_DP10.vcf > bathc1_75_AD_SNPs_GQ50_AD08_DP10_D.vcf
grep '_D' batch1_75_AD_SNPs_GQ50_AD08_DP10.vcf | grep -v '#' >> bathc1_75_AD_SNPs_GQ50_AD08_DP10_D.vcf

d='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI/batch1_86_D'
cd $d
bcftools view -h filtered_SNPs.filtered.vcf > batch1_86_D_SNPs_GO50_AD08_DP10_D.vcf
bcftools view -H filtered_SNPs.filtered.vcf | grep '_D' >> batch1_86_D_SNPs_GO50_AD08_DP10_D.vcf
