#!/bin/bash
cp ./intermediatedata/BRalignment2bp.afa ./intermediatedata/BRalignment2bp.afa.bkp
mafft --auto --thread -1 --keeplength --addfragments ./intermediatedata/BRalignment2bp.afa ./inputdata/NewSeqs.fasta > ./intermediatedata/BRalignment2bp.afa
