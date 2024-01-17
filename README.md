# Germline pipeline V7 Dec 2023

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Feature files](#feature_files)
3. [Get feature files](#get_feature_files)
4. [Detect positions with candidate GVs based on pools samples](#pools_analysis)
5. [Detect positions with GVs based on individual samples](#samples_analysis)
6. [Statistical tests data](#statistical_tests)
7. [Statistical tests filters](#statistical_filters)
8. [Combine filtered files](#combine_filtered_files)
9. [Adding QC indications to variants file](#qc_indications)
10. [Conversion from in-house variants format to standard vcf file](#vcf_conversion)
11. [Adding biological annotations](#bcftools_annotations)
12. [Filtering vcf files using INFO tag data with bcftools](#bcftools_filters)
13. [Intersection with repeats db and Sane genome regions](#bedtools_intersection)
14. [Variant clustring based on location and gnomAD frequency](#variant_clustering)
15. [Filtering only cluster representatives as the filnal list](#cluster_representatives_extraction)
16. [Converting to DS list format](#ds_format)
17. [Collect data and stat from the DS tsv file](#ds_data_extraction)

This repository include scripts and instructutions to create germline analyses based on Nucleix NGS features for individual libraries and pools.
The pipeline scripts are expected to work in AWS linux envirounment with Slurm HPC system configued.

## prerequisites <a name="prerequisites"></a>
Before the Python scripts below can be applied, several prerequisites should be fullfilled

### Python virtual environment should exist with the following non-standard libaries:

The following packages should be inluded in the virtual envirounment
```
biopython==1.81
importlib-metadata==6.8.0
memory-profiler==0.61.0
numpy==1.26.0
pandas==2.1.1
psutil==5.9.5
python-dateutil==2.8.2
pysam==0.21.0
pytz==2023.3.post1
scipy==1.11.3
six==1.16.0
tzdata==2023.3
visidata==2.11.1
zipp==3.17.0
```
The envirounment should be activated before running the scripts below

```
source activate [virtual env name]/bin/activate
```

### Feature files <a name="feature_files"></a>

Consensus Nucs, Consensus Ins, Consensus Del for all samples and pools used. All feature files should be opened (unzipped) and located within parent features directory with subdirectory for each sample:

```
features
├── Exp101G.VU_LNG022.plasma.digested
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Del
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Ins
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Nucs
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr2.Consensus.Del
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr2.Consensus.Ins
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr2.Consensus.Nucs
:
├── Exp101G.VU_LNG022.plasma.digested
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Del
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Ins
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Nucs
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Del
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Ins
│   ├── Exp101G.VU_LNG022.plasma.digested.Chr1.Consensus.Nucs
:

```
## Slurm script to get feature files from the s3 bucket and unzip them <a name="get_feature_files"></a>

In AWS with slurm installed, the follwing script create the above tree structure of feature files. Features are copied from their location on the nucleix.features bucket, unzipped and placed under their sample (library) subdir in the local directory (aws EFS).

```
python slurm_transfer_files.py -s [samples list] -fd [features directory]
```

Example format for samples list:

```
WGS001N.SL_CON666.plasma.digested
WGS001N.SL_CON833.plasma.digested
WGS001N.SL_CON1188.plasma.digested
WGS001N.SL_CON1023.plasma.digested
WGS001N.SL_CON1203.plasma.digested
WGS001N.SL_CON471.plasma.digested
WGS001N.SL_CON334.plasma.digested
WGS001N.SL_CON662.plasma.digested
WGS001N.SL_CON884.plasma.digested
WGS001N.SL_CON874.plasma.digested
:

```

An example command from the v7 run in Dec 2023:

```
python slurm_transfer_files.py -s case_samples_v7.txt -fd /data/users/erane/germline/features/
```


## Detect positions with suspected variability based on pools samples <a name="pools_analysis"></a>

The following script uses the slurm grid to create lists of positions with suspected variability, based on pool samples. The lists of pools samples should include in the second column, the number of samples included in the pool (see below)
```
python slurm_variant_pools.py -outdir [output directory] -c [control_pool_samples_file] -t [case_pool_samples_file]
```


In version 7 analysis (Dec 2023), the following command was applied:

```
python slurm_variant_pools.py -outdir /data/users/erane/germline/variants_pools_v7
```

The pool lists (with number of samples) used for v7 runs. They are identical to the v4 pools:
control_pool_samples_v7.txt:
```
Exp241N.normalPool.WB.untreated	40
Exp900V.controls.plasma.digested	102
Exp65G.normalPool2.plasma.digested	28
Exp31G.normalPool8.plasma.digested	11
Exp31G.normalPool7.plasma.digested	10
Exp31G.normalPool6.plasma.digested	19
Exp31G.normalPool5.plasma.digested	14
Exp31G.normalPool4.plasma.digested	10
Exp31G.normalPool3.plasma.digested	9
Exp31G.normalPool2.plasma.digested	11
Exp31G.normalPool1.plasma.digested	14
Exp7W.normalPool.WB.untreated	10
CMP003N.LOD_050623_4.plasma.untreated	40
WGS900V.controls.plasma.digested	158
TESExp900V.controls.plasma.digested	44
```
case_pool_samples_v7.txt:
```
Exp900V.cancer.plasma.digested	105
Exp230N.normalPool.NT.untreated	40
TESExp900V.cancer.plasma.digested	39
WGS900V.cancer.plasma.digested	108

```
This script was changed in this version to make also the sorting. The output files "*_sorted.tsv" are therefore already sorted and there is no need for a subsequent sorting script.
Another conceptual difference in the detected postion is that the positions set includes also positions for which at least one sample is non-reference homozygous, in order to get in the final positions set also positions which include only non-reference genotypes across all samples.

Overall number of position with suspected variablity in V7 analysis: 

123685010

## Detect positions with variability based on individual case and control samples <a name="samples_analysis"></a>

The following script uses the slurm grid to add to the lists observed variability in individual case and control samples. See lists of samples which were used in the v7 analysis in this directory.
Overall in version 7 there were 484 control samples and 261 case samples.
This stage doe snot filter anything just add information to the previous candidate positions detected in the pools stage.

```
python slurm_detect_variants_v7.py -outdir [output directory] -c [control_pool_samples_file] -t [case_pool_samples_file]
```
Note the slurm script call the following script in the grid:

```
detect_variants_v7.py
```
in parallel fashion for each chunk of bases.
The case and control samples are hard coded in the arguments section of detect_variants_v7.py and are currently:

```
control_samples_v7.txt
case_samples_v7.txt
```

In version 7 analysis (Dec 2023), the following command was applied:

```
python slurm_detect_variants_v7.py -outdir /data/users/erane/germline/variants_genotypes_v7 -posdir /data/users/erane/germline/variants_pools_v7
```

The total number of positions with variability in at least one sample in the analysis:
```
full analysis: 52548542
partial analysis: 43792889
```

## Add statistical test for differences between cases and controls <a name="statistical_tests"></a>

The following script uses the slurm grid to take lists of potential GV detected by the previous step and make statistical tests. Tests include binomial and propotions test between case and controls, binomial test between case and gnomAD external DB. The stat tests p-values are added to output file and filtering is done in the next steps. Also, order test is conducted for the variant frequencies (case variant allele frequency should not be between the control and the gnomAD frequencies) and coverage test to ensure minimal mean coverage at that position. See exact definitions and threholds for the tests in the germline pipeline document ('germline pipeline 151123.pptx')
This script also change binary genotype (BG) of samples with low coverage to '.|.'.
This step adds statistical information but does not filter based on them. The only applied filter is for minimal overall number of total alternative alleles in all samples (-min_allele_coun argument of filter_variants_freqs_v7.py)

```
python slurm_filter_variants_freqs_v7.py -workdir [workdir]
```

This script calls the following script within the AWS grid nodes:

```
filter_variants_freqs_v7.py
```


The case and control samples are hard coded in the arguments section of filter_variants_freqs_v6.py and are currently:
```
control_samples_v7.txt
case_samples_v7.txt
```

In version 7 analysis (Dec 2023), the following command was applied:
```
python slurm_filter_variants_freqs_v7.py -workdir /data/users/erane/germline/variants_filtered_v7/
```



## Apply statistical filters and mean coverage filters to retain relevant GVs <a name="statistical_filters"></a>

The following script uses the slurm grid to take lists of potential GV with statistical test results added by the previous step and apply filters to retain 'interesting' variants. Tests which applied in this step are: (1) p-value of propotions test between case and controls < E10-3 (2) p-value of binomial test between case and gnomAD3.1 frequency < E10-10 (3) order test: case variant allele frequency cannot be between the control and the gnomAD frequencies and the gnomad frequency should be closer to the controls than the case frequencies in a 1:2 ratio (4) mean coverage test - the mean coverage of a given position across all samples should be higher than 10

```
python slurm_filter_variants_freqs_reduce_v7.py -workdir [workdir]
```

This script call the following script in the AWS grid nodes:
```
filter_variants_freqs_reduce_v7.py 
```

The case, control and pool samples should be provided and are hard coded in the arguments section of filter_variants_freqs_reduce_v7.py and are currently:
```
control_samples_v7.txt
case_samples_v7.txt
control_pool_samples_v7.txt
case_pool_samples_v7.txt
```
For the partial analysis (without WGS samples besides EG controls), the following samples were used 
```
control_samples_v7p.txt
case_samples_v7p.txt

```

In version 7 analysis (Dec 2023), the following command was applied:

```
python slurm_filter_variants_freqs_reduce_v7.py -workdir /data/users/erane/germline/variants_filtered_v7/
```




## Combining files to create a single file of variants which pass the statistical and coverage tests from all chunks <a name="combine_filtered_files"></a>

Within the working directory where the variants which passed the previous step are located, type:
```
head -1 2_0_10000000_stat_reduced_v7.vcf > header.vcf
awk FNR-1 *_stat_reduced_v7.vcf > all_stat_reduced_tmp.vcf
cat header.vcf all_stat_reduced_tmp.vcf > all_stat_reduced.vcf
rm all_stat_reduced_tmp.vcf
```

When running the statistical filters step on the slurm grid slurm log files are created (slurm*.out), which contain standard outputs of the filter_variants_freqs_reduce_v7.py scripts. This output reports the number of variants which fail in each particular test.

The following command can be used to combine these outputs to a single file

```
cat slurm-*out | grep -v chunk | grep -v return | grep -v RuntimeWarning | cat slurm_stat_header.txt - >  all_slurms.txt
```

The following python script report summary for the number of variants which fail in each test in all chunks together based on this file:

```
python pass_count_report.py -i all_slurms.txt
```
The program print report to standard output, in the following format:

```
total variants                                                         25300907
case_population_test_fail                                              23335562
case_control_test_fail                                                 25250749
coverage_control_test_fail                                               357373
coverage_case_test_fail                                                  341558
af_order_test_fail                                                      7882404
pass variants                                                              9336
```


## Additional QC indications <a name="qc_indications"></a>

The following script takes as input the combined output file from the previous step (all_stat_reduced.vcf) and add several QC indications 

```
python qc_indications_v7.py -i [variants file] -case [case samples list with backup samples if available] -control [control samples list with backup samples if available] -pools_case [case pools] -pool_control [control pools] -o [output variants file with QC data] 
```

The following command was used, for example in V7:

```
python qc_indications_v7.py -i all_stat_reduced.vcf -o variants_reduced_v7.vcf -case case_samples_with_backup_v7.txt -control control_samples_v7.txt
```


## Conversion to standard vcf4.3 format <a name="vcf_conversion"></a>

The following scripts takes an output format of the previous stage (could be gzipped or unzipped), no need to be indexed, and convert to standard vcf4.3 format which is useful for manipulations with public tools, includiung filtering with bcftools view and annotations with bcftools csq
Note that the program will output all variants types at secific positions, and therefore might be several output lines for each input line. Usually only one of these is statistically significant and it should be retained in subsequent filtering. A change from the previous versions is that in the header tags whoch will be used during variant clustering are being added
This command was tested on the local machine and not in the AWS. 

```
python convert2vcf_v7.py  -i [input_tsv_variants_list] -o [output_vcf_file_name] -control [list of control samples] -case [list of case samples] -local
```

The following command was used, for example, in V7:
```
python convert2vcf_v7.py -i variants_reduced_v7.vcf  -o variants_reduced_formal_v7.vcf -case case_samples_v7.txt -control control_samples_v7.txt -local
```

The script has an option to create per sample information (format section of formal vcf file), which include number of alt consensus reads, overall consensus reads and genotype. To run the script to create extensive format with information for each genotype. Example run in v7:
```
python convert2vcf_v7.py -i variants_reduced_v7.vcf  -o variants_reduced_formal_v7.vcf -case case_samples_v7.txt -control control_samples_v7.txt -control control_samples_no_eg_v7p.txt  -case case_samples_v7p.txt -all_samples -local -full
```

Note that that it is important to accurately provide the -control and -case files, as this annotation dicaate the calculation of case/control allele frequencies in the vcf
The -all_samples option tells the script to report per-sample information in the format section for all samples, and not only to the samples which appear in the -control and -case files


## Adding biological annotations <a name="bcftools_annotations"></a>

Biological annotation including host genes and effect of variant (coding, amino acid change, etc.) were added by csq program of bcftools (https://samtools.github.io/bcftools/howtos/csq-calling.html)

```
bcftools csq [input vcf file] -f [genome path] -g [annotations path] -p a > [output_annotated.vcf]
```
Note that the -p a indicates that the data is not phased. We don't know if adjucent variants falls on the same chromosome or on the sister chromosome.

The following command was used, for example, in V7:

```
bcftools csq all_new_stat_p_plusinfo_formal_full.vcf.bgz  -f ~/Genomes/hg38.fa -g ~/Genomes/hg38.gff3.gz -p a > all_new_stat_p_plusinfo_formal_full_annotated.vcf
bcftools csq variants_reduced_v7_new_formal.vcf.bgz  -f ~/Genomes/hg38.fa -g ~/Genomes/hg38.gff3.gz -p a > variants_reduced_new_formal_annotated_v7.vcf

```

## Variant filtering command examples with bcftools view <a name="bcftools_filters"></a>

Filters based on variety of QC criteria. 

Differences in filtering comparing to V4:
 <li> The case-control is now based on propotions test (PR tag) and the final thershold set to < 1*10-3</li>
 <li> Variants with no calls in gnomAD are excluded </li>
 <li> Maximum homoplymer stretch is 5 (instead 6 before)</li>
 <li> Variants with Hardy–Weinberg test p-value with p-value (PH tag) < 1*E-3 are excluded</li>
 <li> Variants with no information in gnomAD (or gnomaAD frequency is set to 8*E-16</li>
<br>

```
bcftools view -e 'NF_CONT > 0.02 || NF_CASE > 0.02 || abs(AFRF_CASE - AFRF_CONT) > 0.3 || ((VF_CASE < 0.8) && (VN_CASE > 0)) || ACR > 2 || ACR < 0.5 || RN > 4 || GA == 8E-6 || PR > 1E-3 || HP > 5 || PH < 1E-3' all_new_stat_p_plusinfo_formal_full_annotated.vcf  > all_new_stat_p_plusinfo_formal_full_annotated_filtered.vcf

bcftools view -e 'NF_CONT > 0.02 || NF_CASE > 0.02 || abs(AFRF_CASE - AFRF_CONT) > 0.3 || ((VF_CASE < 0.8) && (VN_CASE > 0)) || ACR > 2 || ACR < 0.5 || RN > 4 || GA == 8E-6 || PR > 1E-3 || HP > 5 || PH < 1E-3'  variants_reduced_formal_annotated_v7p.vcf.bgz > variants_reduced_formal_annotated_filtered_v7p.vcf

bcftools view -e 'NF_CONT > 0.02 || NF_CASE > 0.02 || abs(AFRF_CASE - AFRF_CONT) > 0.3 || ((VF_CASE < 0.8) && (VN_CASE > 0)) || ACR > 2 || ACR < 0.5 || RN > 4 || GA == 8E-6 || PR > 1E-3 || HP > 5 || PH < 1E-3' all_new_stat_p_no_eg_plusinfo_formal_full_annotated.vcf  > all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered.vcf
```


## Intersection with repeats db and Sane genome regions <a name="bedtools_intersection"></a>

The output of the previous step is intersected with simple repeats file to **exclude** overlaps and then with Sane genome bed file to **retain** more confident regions 
```
bedtools intersect -header -v -a [vcf file] -b /home/eraneyal/Genomes/ucsc_RepeatMasker_hg38_nucleix_sorted_simple.bed | bedtools intersect -header -u -a - -b /home/eraneyal/Genomes/Sane_hg38.bed > [vcf file no repeats and sane]
```

Example commands from the v7 analyses:

```
bedtools intersect -header -v -a variants_reduced_new_formal_annotated_filtered_v7.vcf -b /home/eraneyal/Genomes/ucsc_RepeatMasker_hg38_nucleix_sorted_simple.bed | bedtools intersect -header -u -a - -b /home/eraneyal//Genomes/Sane_hg38.bed > variants_reduced_new_formal_annotated_filtered_no_rep_sane_v7.vcf

bedtools intersect -header -v -a all_new_stat_p_plusinfo_formal_full_annotated_filtered.vcf  -b /home/eraneyal/Genomes/ucsc_RepeatMasker_hg38_nucleix_sorted_simple.bed | bedtools intersect -header -u -a - -b /home/eraneyal/Genomes/Sane_hg38.bed > all_new_stat_p_plusinfo_formal_full_annotated_filtered_no_rep_sane_v7p.vcf

bedtools intersect -header -v -a all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered.vcf  -b /home/eraneyal/Genomes/ucsc_RepeatMasker_hg38_nucleix_sorted_simple.bed | bedtools intersect -header -u -a - -b /home/eraneyal/Genomes/Sane_hg38.bed > all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered_no_rep_sane_v7p.vcf
```

## Clustering of variants based on close proximity and gnomAD variant frequency <a name="variant_clustering"></a>

In house script is used for variant clustering. The goal is to reduce effective number of variants for training and QC purposes

```
python cluster_variants.py -i variants_reduced_new_formal_annotated_filtered_no_rep_sane_v7.vcf -o variants_reduced_new_formal_annotated_filtered_no_rep_sane_clustered_v7.vcf
python cluster_variants.py -i variants_reduced_formal_annotated_filtered_no_rep_sane_v7p.vcf -o variants_reduced_formal_annotated_filtered_no_rep_sane_clustered_v7p.vcf
python cluster_variants.py -i lists_debug_Dec27/all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered_no_rep_sane_v7p.vcf -o lists_debug_Dec27/all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered_no_rep_sane_clustered.vcf
```

## Extraction of the cluster representatives as the final list <a name="cluster_representatives_extraction"></a>

```
bcftools view -i 'RP==1' variants_reduced_new_formal_annotated_filtered_no_rep_sane_clustered_v7.vcf  > variants_reduced_new_formal_annotated_filtered_no_rep_sane_cluster_reps_v7.vcf
bcftools view -i 'RP==1' variants_reduced_new_formal_annotated_filtered_no_rep_sane_clustered_v7p.vcf  > variants_reduced_new_formal_annotated_filtered_no_rep_sane_cluster_reps_v7p.vcf
bcftools view -i 'RP==1'  all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered_no_rep_sane_clustered.vcf > all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered_no_rep_sane_clustered_representatives.vcf
```

## Converting to DS list format <a name="ds_format"></a>

The final vcf file is converted to a more friendly tabular format for the DS team. The final list is also arranged such the the frequencies and odds ratios are calculated according to the "minor" allele frequency (taken as the allele frequency in our control samples 

```
python vcf2ds.py -i lists_debug_Dec27/all_new_stat_p_no_eg_plusinfo_formal_full_annotated_filtered_no_rep_sane_clustered_representatives_manual.vcf -o GV_partial_no_eg_analysis_list_Dec23_for_DS.tsv
```
The columns in the output list are:
```
hg38_chr – chromosome according to hg38
hg38_pos – genomic position according to hg38
nucleix_pos – genomic position according to Nucleix features. Might be different from hg38_pos for deletions.
ref – reference base in the position according to hg38
alt - alternative base/indel which might be found in that site
gv – the germline variant, defined as the minor allele genotype in the controls cohort at that position. Could be the alt base (most often) or the reference base.
annotation – annotations based on bcftools csq program. This field is populated only for variants within proteins
p_control_gnomad – p-value of binomial test between the control gv allele frequency and gnomAD frequency.
p_case_gnomad - p-value of binomial test between the cases gv allele frequency and gnomAD frequency.
p_case_control_binomial - p-value of binomial test between the cases gv allele frequency and controls gv allele frequency.
p_case_control_prop - p-value of proportions test between the cases gv allele frequency and controls gv allele frequency.
gv_alleles_control – number of gv alleles in the control subjects of the training cohort
total_alleles_control – total number of alleles in the control subjects of the training cohort (2 x number of control subjects)
af_control – gv allele frequency in the control subjects of the training cohort
gv_alleles_case - number of gv alleles in the case subjects of the training cohort
total_alleles_case - total number of alleles in the case subjects of the training cohort (2 x number of case subjects)
af_case - gv allele frequency in the case subjects of the training cohort
odds_ratio - af_case / af_control
gnomad_v3.1_af - gv frequency in the population database gnomAD v3.1
for each sample: number of gv alleles. Could be 0,1,2 or -1 for no data. -1 should appear in this table only for training samples, so no need to impute values of validation subjects.
```



## Collect data and stat from the DS tsv file <a name="ds_data_extraction"></a>

The following scripts were used to extract information from teh DS tsv file

The following example collects odds ratios and gnomad frequencies from the table and count SNVs, insertions and deletions
```
python ds_list_stat.py -i GV_partial_analysis_list_Dec23_for_DS.tsv > partial_gnomad_odd_ratios.tsv
```

The following example collects from the list odds ratios in the training samples and another test subset

```
python check_wgs_ratios_ds_format.py -i GV_partial_analysis_list_Dec23_for_DS.tsv -o train_test_odds_ratios.txt
```

The last example produces a file with values of the odds_ratio in the training set and test set. The following awk command allow you to get the total number of lines in which there is agreement in the direction of the tarin and test odds_artios (both smaller than one or larger than one): 

```
awk 'BEGIN {sum=0} {if (log($1) * log($2) > 1) sum+=1} END {print(sum)}' train_test_odds_ratios.txt
```





