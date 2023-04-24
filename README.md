# LRNRisk
Materials for the LRNRisk project

## lrnrisk_prototype.py
Python script producing LRNRisk HTML file for a single genome
```
usage: lrnrisk_prototype.py -v </path/to/vfdb_output.tab> -a </path/to/pima_amr_output.tab> -o </path/to/output/directory/> [other options]

optional arguments:
  -h, --help            show this help message and exit
  -g GTDB, --gtdb GTDB  Path to gtdbtk.bac120.summary.tsv, the tab-separated file containing GTDB-Tk results
  -v VIRULENCE, --virulence VIRULENCE
                        Path to tab-separated file containing VFDB virulence factors detected via BLAST
  -a AMR, --amr AMR     Path to tab-separated file containing Pima AMR determinants detected via BLAST
  -b BLACKLIST, --blacklist BLACKLIST
                        Path to tab-separated file containing blacklisted high-risk virulence factors
  -dv VF_DISTRIBUTION, --vf_distribution VF_DISTRIBUTION
                        Path to tab-separated file containing virulence factor distribution
  -da AMR_DISTRIBUTION, --amr_distribution AMR_DISTRIBUTION
                        Path to tab-separated file containing AMR determinant distribution
  -o OUTPUT, --output OUTPUT
                        Path to desired output HTML file
```

## lrnrisk_test.sh
Sample command to run lrnrisk_prototype.py:
```
python lrnrisk_prototype.py --gtdb gtdb_example/gtdbtk.bac120.summary.tsv --virulence pima_example/pima_output/vfdb_core.fasta.tabular --amr pima_example/pima_output/AMR_pima_md_2023_02_02.fasta.tabular --blacklist lrnrisk_prototype_blacklist.tsv --vf_distribution vfdb_dist_beta_80p80q.tsv --amr_distribution resfinder_dist_beta_80p80q.tsv --output test.html
```

## lrnrisk_logo.jpg
LRNRisk logo (displayed when users open HTML output file)

## lrnrisk_prototype_blacklist.tsv
"Blacklisted genes"; used to flag a user-supplied genome as "high-risk" if one or more genes is detected in the query genome; supplied as input to `lrnrisk_prototype.py --blacklist`

## vfdb_dist_beta_80p80q.tsv
Virulence factors detected in publicly available genomes, broken down by species; supplied to `lrnrisk_prototype.py --vf_distribution` 

## resfinder_dist_beta_80p80q.tsv
Antimicrobial resistance (AMR) determinants detected in publicly available genomes, broken down by species; supplied to `lrnrisk_prototype.py --amr_distribution`

## gtdb_example/gtdbtk.bac120.summary.tsv
Genome Taxonomy Database (GTDB) Toolkit (GTDB-Tk) output file for a sample genome; supplied to `lrnrisk_prototype.py --gtdb`

## pima_example/pima_output/vfdb_core.fasta.tabular
Virulence Factor Database (VFDB) BLAST results for a single genome; supplied as input to `lrnrisk_prototype.py --virulence`

## pima_example/pima_output/AMR_pima_md_2023_02_02.fasta.tabular
AMR determinant detection results for a single genome; supplied as input to `lrnrisk_prototype.py --amr`



