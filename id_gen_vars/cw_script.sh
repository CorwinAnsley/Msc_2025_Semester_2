#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l cput=1:00:00
#PBS -l nodes=1:centos7:ppn=1
#PBS -m abe 
#PBS -q bioinf-stud
#PBS -M 2266643a@student.gla.ac.uk

plink --bfile /export/biostuds/2266643a/BIOL5300_genevarid/cw_data/2266643/quantitative/2266643_2 \
--linear \
--maf 0.01 \
--hwe 0.001 \
--mind 0.05 \
--geno 0.99 \
--out /export/biostuds/2266643a/BIOL5300_genevarid/cw_results/cw_output_qc__chr_2
