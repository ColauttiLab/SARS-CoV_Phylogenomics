#!/bin/bash

date

module load gcc/8.3.0
module load r/4.0.0

Rscript --vanilla Distcalc.R

date
