
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
cloud_vcf='CryptoD_test.genotype.filter.haploid.cloud.vcf'
merged_vcf='CryptoD_test.genotype.filter.haploid.merged.cloud.vcf'
gsutil cp gs://4b66fc8a-tmp/cromwell-executions/gatk_process_cohort/f568dfdf-4c93-4c3b-b93d-c48cc0680aca/call-gatk_filter_genotypes_task/CryptoD_test.genotype.filtered.vcf $cloud_vcf
java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar -T SelectVariants -env \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta \
  -V $cloud_vcf -o cloud_SNPs.vcf -se '.+variant$'

java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar -T SelectVariants -env \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta \
  -V $cloud_vcf -o cloud_INDELs.vcf -se '.+variant2$'

# Process cloud
java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar -T CombineVariants \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta \
  --assumeIdenticalSamples --variant cloud_SNPs.vcf --variant cloud_INDELs.vcf \
  -o $merged_vcf -genotypeMergeOptions UNSORTED

# calculate concordance
java -jar /xchip/gtex/xiaoli/tools/GenomeAnalysisTK.jar \
  --comp CryptoD_test.genotype.filter.onprem.vcf.gz \
  --eval $merged_vcf \
  -R /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta \
  -T GenotypeConcordance  -o genotConcord.txt

#
bcftools query -f '%CHROM\_%POS\_%REF\_%ALT{0} %AF %AC %AN\n' $merged_vcf > cloud_AF.txt
bcftools query -f '%CHROM\_%POS\_%REF\_%ALT{0} %AF %AC %AN\n' CryptoD_test.genotype.filter.onprem.vcf.gz > onprem_AF.txt
