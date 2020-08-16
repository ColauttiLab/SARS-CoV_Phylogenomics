# COVID Analysis Pipeline

Analysis Files from 2020 COV-19 genome study by Sjaarda et al:
"Chasing the origin of SARS-CoV-2 in Canada’s COVID-19 cases: A genomics study"
https://www.biorxiv.org/content/10.1101/2020.06.25.171744v1

# scripts Overview

1. Download latest nextstrain squence data from GISAID (https://www.gisaid.org)
  * OR to reproduce the published analysis, download the archived sequence data from <<LINK>>
2. Run makeBlastDB.sh to create a local fasta database for the nextseq sequences (this will take a while)
3. Run NextStrainSetup.R to BLAST reference files to nextstrain sequences and subset analysis
4. Run mafftalign.sh to create an alignment of the full set of sequences

## Instructions
1. You may need to change your local path in scripts/retrieveSeq.sh to point to blast (especially if you aren't using a Mac) 

## TO DO:

- [ ] Add all Wuhan strains
- [ ] Use updated Wuhan root: Wuhan/WH04/2020
- [ ] check/clean up libraries in NextStrainSetup.R & NextStrainAnalysis_Main.R
- [X] retrieveSeq.sh -- replace directories with generic names
- [ ] Add code for phylogeny
  - [ ] Add pangolin lineage names to phylogeny (A, B, B1, B1.5, etc.)
- [ ] Add link to data archive
- [ ] Add citation to paper
  
