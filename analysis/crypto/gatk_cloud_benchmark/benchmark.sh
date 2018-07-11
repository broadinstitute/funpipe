
# onprem merge
cp call-FilterIndels/execution/filtered_INDELS.vcf* /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI/pilot_4_samples_test/
cp call-FilterSnps/execution/filtered_SNPs.vcf* /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI/pilot_4_samples_test/

java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar \
  --variant filtered_INDELS_D.vcf --variant filtered_SNPs_D.vcf \
  -o CryptoD_test.genotype.filter.onprem.vcf -genotypeMergeOptions UNSORTED \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta
  -T CombineVariants --assumeIdenticalSamples

# cloud data
# separate
java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar -T SelectVariants -env \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta \
  -V /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/test_cloud/output/CryptoD_test.genotype.filter.0608.vcf \
  -o cloud_SNPs.vcf \
  -se '.+variant$'

java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar -T SelectVariants -env \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta \
  -V /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/test_cloud/output/CryptoD_test.genotype.filter.0608.vcf \
  -o cloud_INDELs.vcf \
  -se '.+variant2$'

# Process cloud
java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar -T CombineVariants -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta --assumeIdenticalSamples --variant cloud_SNPs.vcf --variant cloud_INDELs.vcf -o CryptoD_test.genotype.filter.cloud.vcf -genotypeMergeOptions UNSORTED

# calculate concordance
java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar \
  --comp CryptoD_test.genotype.filter.onprem.vcf.gz \
  --eval CryptoD_test.genotype.filter.cloud.vcf \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta \
  -T GenotypeConcordance  -o genotConcord.txt
