#!/bin/bash
#$ -cwd           
#$ -pe smp 1      
#$ -l h_vmem=16G
#$ -l h_rt=1:0:0  
#module unload java/11.0.2
module load java/17.0.0

java -Djava.io.tmpdir=/data/.../Exomiser_Scripts -jar /data/WHRI-Phenogenomics/projects/pheval/configurations/exomiser-14.0.0/2406/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar \
--analysis /data/home/hhz220/ndd_exomiser/preset-phenotype-analysis.yml \
--output-filename NDD_HiPhive_Prioritiser \
--output-format TSV_GENE \
--output-directory /data/home/hhz220/ndd_exomiser/output \
--spring.config.location=/data/WHRI-Phenogenomics/projects/pheval/configurations/exomiser-14.0.0/2406/application.properties
