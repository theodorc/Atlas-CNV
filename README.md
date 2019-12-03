# Atlas-CNV

## Overview
Atlas-CNV is a method to identify copy number variants (CNVs) using Read Depth coverage at the exon level in the next-generation sequencing (NGS) of Gene Targeted Panels. Target Exons of CNVs are detected within a batch of samples (47) from the same probe capture experiment (midpool). Individually detected exons are C-scored (Z-like scores) which allows for prioritization of high quality CNVs especially for single-exon deletions or duplications. Gene and exon bar plots are also produced for additional evaluation.

Please cite: 
Chiang, T., Liu, X., Wu, T. et al. Atlas-CNV: a validated approach to call single-exon CNVs in the eMERGESeq gene panel. Genet Med 21, 2135–2144 (2019) doi:10.1038/s41436-019-0475-4

## Dependencies
Initial release used 

* Perl v5.12.2 
   * use strict;
   * use Getopt::Long;
   * use Pod::Usage;
   * use Config::General;

* R v3.1.1. 
   * library(reshape2)
   * library("optparse")
   * library(plotrix)

## Usage
```
./atlas_cnv.pl
Usage:
    atlas_cnv.pl [options] [file ...]

    Options:

     --config             path to the configuration file [default = config]
     --panel              path to the panel design file (4 tab-delimited columns)
                          Exon_Target       Gene_Exon       Call_CNV        RefSeq
                          1:1220087-1220186 SNP_1           N               rs2144440

     --sample             path to the CNV sample tab-delimited file containing sample_id, gender, midpool
     --threshold_del      log2 ratio threshold to call a deletion [uses defaults in atlas_cnv.R]
     --threshold_dup      log2 ratio threshold to call a duplication [uses default in atlas_cnv.R]
     --help|?             prints a brief help message
     --man                prints a man page

Options:
    --config-path
            The path to the file containing configuration items.

    --panel The path to the panel design tab-delimited file. Exon_Target,
            Gene_Exon, Call_CNV, RefSeq

    --sample
            The path to the tab-delimited file containing sample_id, gender,
            and midpool

    --threshold_del
            The log2 ratio threshold to call a deletion log2(sample/median)
            [uses defaults in atlas_cnv.R]

    --threshold_dup
            The log2 ratio threshold to call a duplication
            log2(sample/median) [uses defaults in atlas_cnv.R]

    --help  Print a brief help message and exits.

    --man   Prints the manual page and exits.

```

## Example: *sample* and *panel* file format
panel file (exact header required):
```
Exon_Target          Gene_Exon      Call_CNV  RefSeq
1:1220087-1220186    SNP_1          N         rs2144440
1:3083663-3083762    SNP_2          N         rs2651899
1:3611843-3611942    SNP_3          N         rs3765731
1:6279321-6279420    RNF207-001_18  N         rs846111
1:8487274-8487373    SNP_4          N         rs301797
1:11850737-11850955  MTHFR-001_11   Y         NM_005957_cds_0
1:11851264-11851383  MTHFR-001_10   Y         NM_005957_cds_1
1:11852335-11852436  MTHFR-001_9    Y         NM_005957_cds_2
1:11853964-11854146  MTHFR-001_8    Y         NM_005957_cds_3
etc...
```
sample file (no header): *sample*  *sex*  *midpool*
```
1000011111_HFG5NBCXY-1-IDMB11  F  IRC_MB-MID00209.CLEMRG_1pA
1000022222_HFG5NBCXY-1-IDMB22  M  IRC_MB-MID00209.CLEMRG_1pA
1000033333_HFG5NBCXY-1-IDMB33  F  IRC_MB-MID00210.CLEMRG_1pA
1000044444_HFG5NBCXY-1-IDMB44  F  IRC_MB-MID00210.CLEMRG_1pA
1000055555_HFG5NBCXY-1-IDMB55  M  IRC_MB-MID00209.CLEMRG_1pA
etc...

```
