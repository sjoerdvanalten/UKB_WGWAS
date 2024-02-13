detailed lookup

#wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/categorical-20002-both_sexes-1222.tsv.bgz
#
#zcat categorical-20002-both_sexes-1222.tsv.bgz | grep 83998003 #similar negative effect size in pan-ukb, p-val of 0.35
#zcat categorical-20002-both_sexes-1222.tsv.bgz | grep 154123670
#zcat categorical-20002-both_sexes-1222.tsv.bgz | grep 83998003 #similar negative effect size in pan-ukb, p-val of 0.35
#zcat categorical-20002-both_sexes-1222.tsv.bgz | grep 83998003 #similar negative effect size in pan-ukb, p-val of 0.35

#wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000181/GCST90000181_buildGRCh38.tsv

wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010681/GCST010681_buildGRCh37.tsv.gz

zcat GCST010681_buildGRCh37.tsv.gz | grep 83998003
zcat GCST010681_buildGRCh37.tsv.gz | grep 154123670
zcat GCST010681_buildGRCh37.tsv.gz | grep 20834175

#breast cancer 
wget https://bcac.ccge.medschl.cam.ac.uk/files/oncoarray_bcac_public_release_oct17.txt.gz #file retrieved from https://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-results-breast-cancer-risk-2017
zcat oncoarray_bcac_public_release_oct17.txt.gz | grep rs2306412


The breast cancer genome-wide association analyses were supported by the Government of Canada through Genome Canada and the Canadian Institutes of Health Research, 
the ‘Ministère de l’Économie, de la Science et de l’Innovation du Québec’ through Genome Québec and grant PSR-SIIRI-701, The National Institutes of Health (U19 CA148065, X01HG007492), 
Cancer Research UK (C1287/A10118, C1287/A16563, C1287/A10710) and The European Union (HEALTH-F2-2009-223175 and H2020 633784 and 634935). All studies and funders are listed in Michailidou 
et al (Nature, 2017).

#check in Ben-Neales:
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/E10.gwas.imputed_v3.both_sexes.tsv.bgz -O E10.gwas.imputed_v3.both_sexes.tsv.bgz #T1D
zcat E10.gwas.imputed_v3.both_sexes.tsv.bgz | grep 154123670
zcat E10.gwas.imputed_v3.both_sexes.tsv.bgz | grep 20834175
zcat E10.gwas.imputed_v3.both_sexes.tsv.bgz | grep 83998003


wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/C50.gwas.imputed_v3.female.tsv.bgz -O C50.gwas.imputed_v3.female.tsv.bgz #Breast Cancer 

