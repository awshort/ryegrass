#List of commands used to modify the VCF and perform analysis on LG5 of the ryegrass genome
#Filter phased VCF to just LG5. Remove indels. Filter for biallelic sites. Minor allele freqeuncy cutoff of 0.005
vcftools --gzvcf phased.output_merged_new_ID_miss075_dp10.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.005 --recode --out no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10
#bgzip VCF
bgzip no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#index VCF
tabix no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz
#Filter to the 7 LGs
bcftools view 'no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' --regions 1_RagTag,2_RagTag,3_RagTag,4_RagTag,5_RagTag,6_RagTag,7_RagTag -o 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Calculate Fst between susceptible and resistant samples
vcftools --vcf 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf \
--weir-fst-pop s_pops.csv \
--weir-fst-pop r_pops.csv \
--out 1_7_RagTag.no_indels.maf_0.005.biallelic.fst.ryegrass_s_r 
#VCF with just resistant samples
vcftools --keep r_pops.csv --vcf 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf --recode --out r.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10
#VCF of just susceptible samples
vcftools --keep s_pops.csv --vcf 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf --recode --out s.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10
#calculate Fst between pops
vcftools --vcf 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf \
--weir-fst-pop SLB_pops.csv \
--weir-fst-pop PR_pops.csv \
--out 1_7_RagTag.no_indels.maf_0.005.biallelic.fst.ryegrass_SLB_PR 

#calculate Fst between pops
vcftools --vcf 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf \
--weir-fst-pop GULF.csv \
--weir-fst-pop L60.csv \
--out 1_7_RagTag.no_indels.maf_0.005.biallelic.fst.ryegrass_GULF_L60

#bgzip VCF file
bgzip r.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf

#Index VCF file
tabix r.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz

#Filter to chromsome 5
bcftools view 'r.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' --regions 5_RagTag -o r.5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Bgzip VCF file 
bgzip s.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Index VCF file
tabix s.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz
#Filter to chromosome 5
bcftools view 's.1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' --regions 5_RagTag -o s.5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf

#Set genomics general script to the path
export PYTHONPATH=$PYTHONPATH:'/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master'
#Activate conda environment with python 2
source activate py2 
#Convert VCF to geno file
python '/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master/VCF_processing/parseVCF.py' -i '1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf' --skipIndels > 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.geno
#calculate Fst between susceptible and resistant in windows
python2.7 '/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master/popgenWindows.py' -g 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.geno -w 100000 -m 100 -s 1000 -f phased -p S -p R --popsFile pop_pixy_2.txt -T 10 -o Fst.w_10000.m100.s_1000.ryegrass.csv --writeFailedWindow
#Set genomics general script to the path
export PYTHONPATH=$PYTHONPATH:'/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master'
#Activate conda environment with python 2
source activate py2 
#calculate Fst between susceptible and resistant in windows
python2.7 '/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master/popgenWindows.py' -g 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.geno -w 100000 -m 100 -f phased -p S -p R --popsFile pop_pixy_2.txt -T 10 -o Fst.w_100000.m100.ryegrass.csv --writeFailedWindow
#calculate Fst between susceptible and resistant in windows
python2.7 '/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master/popgenWindows.py' -g 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.geno -w 100000 -m 100 -f phased -p GULF -p L31 -p L46 -p L60 -p PR -p SLB --popsFile pop_pixy_3.csv -T 10 -o Fst.w_100000.m100.ryegrass_pops.csv --writeFailedWindow
#Set genomics general script to the path
export PYTHONPATH=$PYTHONPATH:'/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master'
#Activate conda environment with python 2
source activate py2 
#calculate Fst between susceptible and resistant in windows
python2.7 '/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master/popgenWindows.py' -g 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.geno -w 25000 -m 25 -f phased -p S -p R --popsFile pop_pixy_2.txt -T 10 -o Fst.w_25000.m25.ryegrass.csv --writeFailedWindow
#calculate Fst between susceptible and resistant in windows
python2.7 '/media/aidanwilliamshort/extradrive1/Spring_Rotation/genomics_general-master/popgenWindows.py' -g 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.geno -w 25000 -m 25 -f phased -p GULF -p L31 -p L46 -p L60 -p PR -p SLB --popsFile pop_pixy_3.csv -T 10 -o Fst.w_25000.m25.ryegrass_pops.csv --writeFailedWindow

#Bgzip VCF
bgzip 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Index VCF
tabix 1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz
#Output chromosome 5
bcftools view '1_7_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' --regions 5_RagTag -o 5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Bgzip chromosome 5 VCF
bgzip 5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf
#Index chromosome 5 VCF
tabix 5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz
#Create VCF of the highly differentiated region on chromosome 5
bedtools intersect -a '5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf.gz' -b 'high_differentiated_region_sites_coords.bed' -header > high_differentiated_region.5_RagTag.no_indels.maf_0.005.biallelic.phased.output_merged_new_ID_miss075_dp10.recode.vcf

