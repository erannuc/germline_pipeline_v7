# Germline pipeline V7 Dec 2023

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Feature files](#feature_files)
3. [Get feature files](#get_feature_files)
4. [Detect positions with candidate GVs based on pools sample](#pools_analysis)
5. [Detect positions with GVs based on individual samples](#samples_analysis)
6. [Statistical tests data](#statistical_tests)
7. [Statistical tests filters](#statistical_filters)
8. [Combine filtered files](#combine_filtered_files)
9. [Adding QC indications to variants file](#qc_indications)
10. [Conversion from in-house variants format to standard vcf file](#vcf_conversion)
11. [Filtering vcf files using INFO tag data with bcftools](#bcftools_filters)
12. [Adding biological annotations](#bcftools_annotations)

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

Consensus Nucs, Consensus Ins, Consensus Del for all samples and pools used. All feature files should be opened (unzipped) and located within parent features directory with subsirectory for each sample:

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

In AWS with slurm installed, the follwing script create the above tree structure of feature files. Features are copied from their location on the nucleix.features bucket and unzipped and located under their sample (library) subdir in the local directory (aws EFS).

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

An example command from the v7 run in DEc 2023:

```
python slurm_transfer_files.py -s case_samples_v7.txt -fd /data/users/erane/germline/features/
```


## Detect positions with suspected variability based on pools samples <a name="pools_analysis"></a>

The following script uses the slurm grid to create lists of positions with suspected variability, based on pool samples. The lists of pools samples should include in the second column, the number of samples included in the pool (see below)
```
python slurm_variant_pools.py -outdir [output directory] -c [control_pool_samples_file] -t [case_pool_samples_file]
```


In version 4 analysis (Oct 2023), the following command was applied:

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
The case and control samples are hard coded in the arguments section of detect_variants_v4.py and are currently:

```
control_samples_v7.txt
case_samples_v7.txt
```

In version 7 analysis (Dec 2023), the following command was applied:

```
python slurm_detect_variants_v7.py -outdir /data/users/erane/germline/variants_genotypes_v7 -posdir /data/users/erane/germline/variants_pools_v7
```

## Add statistical test for differences between cases and controls <a name="statistical_tests"></a>

The following script uses the slurm grid to take lists of potential GV detected by the previous step and make statistical tests. Tests include binomial test between case and controls, binomial test between case and gnomAD external DB. Note that no filter is conducted (apart of minimal alleles represented in the samples). The stat tests p-values are added to output file and filtering is done in the next step. Also, order test is conducted for the variant frequencies (case variant allele frequency should not be between the control and the gnomAD frequencies) and coverage test to ensure minimal mean coverage at that position. See exact definitions and threholds for lasl tests in the germline pipeline document ('germline pipeline 151123.pptx')
This script also change binary genotype (BG) of samples with low coverage to '.|.'. Low coverage  

```
python slurm_filter_variants_freqs_v4.py -workdir [workdir]
```

This script calls the following script within the grid nodes:

```
filter_variants_freqs_v6.py
```


The case and control samples are hard coded in the arguments section of filter_variants_freqs_v6.py and are currently:
```
control_samples_v4.txt
case_samples_v4.txt
```

In version 4 analysis (Oct 2023), the following command was applied:
```
python slurm_filter_variants_freqs_v4.py -workdir /data/users/erane/germline/variants_filtered_v4_1_1/
```


## Apply statistical filters and mean coverage filters to retain relevant GVs <a name="statistical_filters"></a>

The following script uses the slurm grid to take lists of potential GV with statistical test results added by the previous step and apply filters to retain 'interesting' variants. Tests which applied in this step are: (1) p-value of binomial test between case and < E10-5 (2) p-value of binomial test between case and gnomAD3.1 frequency < E10-13 (3) order test: case variant allele frequency cannot be between the control and the gnomAD frequencies (4) mean coverage test - the mean coverage of a given position across all samples should be higher than 

```
python slurm_filter_variants_freqs_reduce_v4.py -workdir [workdir]
```
The case, control and pool samples should be provided and are hard coded in the arguments section of filter_variants_freqs_reduce_v4.py and are currently:
```
control_samples_v4.txt
case_samples_v4.txt
control_pool_samples_v4.txt
case_pool_samples_v4.txt
```
In version 4 analysis (Oct 2023), the following command was applied:

```
python slurm_filter_variants_freqs_reduce_v4.py -workdir /data/users/erane/germline/variants_filtered_v4_1_1/
```


## Combining files to create a single file of variants which passe the statistical and coverage tests from all chunks <a name="combine_filtered_files"></a>

Within the working directory where the variants which passed the previous step are located, type:
```
head -1 2_0_10000000_stat_reduced_v4.vcf > header.vcf
awk FNR-1 *_stat_reduced_v4.vcf > all_stat_reduced_tmp.vcf
cat header.vcf all_stat_reduced_tmp.vcf > all_stat_reduced.vcf
rm all_stat_reduced_tmp.vcf
```

