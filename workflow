#List of commands used to modify the VCF and perform analysis on LG5 of the ryegrass genome
#Filter phased VCF to just LG5. Remove indels. Filter for biallelic sites. Minor allele freqeuncy cutoff of 0.005
vcftools --gzvcf phased.output_merged_new_ID_miss075_dp10.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.005 --recode --out no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10
#bgzip VCF
bgzip no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#index VCF
tabix no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz
#Filter to the 7 LGs
bcftools view 'no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' --regions 1_RagTag,2_RagTag,3_RagTag,4_RagTag,5_RagTag,6_RagTag,7_RagTag -o 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Filter for MAF of 0.05
vcftools --gzvcf 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz --maf 0.05 --recode --out 1_7_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10
#bgzip VCF
bgzip 1_7_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#index VCF
tabix 1_7_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz
#Filter to LG5
bcftools view '1_7_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' --regions 5_RagTag -o 5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Calculate Fst
vcftools --vcf 5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf \
--weir-fst-pop s_pops.csv \
--weir-fst-pop r_pops.csv \
--out 5_RagTag.no_indels.maf_0.05.biallelic.phased.fst.ryegrass_s_r 
#File with just resistant samples
vcftools --keep r_pops.csv --vcf '5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf' --recode --out r.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10
#File with just susceptible samples
vcftools --keep s_pops.csv --vcf '5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf' --recode --out s.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10
#bgzip VCF
bgzip 5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#index VCF
tabix 5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz
#Filter VCF to highly differentiated region
bedtools intersect -a '5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' -b 'high_differentiated_region_sites_coords.bed' -header > high_differentiated_region.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Create fasta from VCF of the highly differentiated region
'/media/aidanwilliamshort/extradrive2/vcfx_2.0.6b/source/build/vcfx' fasta input='high_differentiated_region.5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf' output=high_differentiated_region.5_RagTag.fasta start=160175000 end=160186000  reference='5_RagTag.LOLMU.fa'
#Create VCF with just the CBP gene 
bedtools intersect -a '5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' -b 'CBP_sites_coords.bed' -header > CBP.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Copy CBP VCF as text file
cp CBP.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf CBP.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.txt
#Use bash script to create header list with A and B haplotypes for each sample
bash '/media/aidanwilliamshort/barbatus1/ryegrass/forloop.bash'
#Create VCF with different coordinates for the highly differentiated region
bedtools intersect -a '5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' -b 'high_differentiated_region_2.bed' -header > high_differentiation_2.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Copy VCF as text file
cp high_differentiation_2.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf high_differentiation_2.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.txt
#Change VCF to FASTA file 
'/media/aidanwilliamshort/extradrive2/vcfx_2.0.6b/source/build/vcfx' fasta input='high_differentiated_region.5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf' output=high_differentiated_region_160175000_160409000.5_RagTag.fasta start=160175000 end=160409000  reference='5_RagTag.LOLMU.fa'
#Create VCF just for defined coordinates in bed file
bedtools intersect -a '5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' -b '5_RagTag_99780000_100070000.bed' -header > 5_RagTag_99780000_100070000.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Change VCF to fasta file
'/media/aidanwilliamshort/extradrive2/vcfx_2.0.6b/source/build/vcfx' fasta input='5_RagTag_99780000_100070000.5_RagTag.no_indels.maf_0.05.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf' output=5_RagTag_99780000_100070000.5_RagTag.fasta start=99780000 end=100070000 reference='5_RagTag.LOLMU.fa'





