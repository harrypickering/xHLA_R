#!/bin/bash
#usage : ./000_master_control_script.sh_id
read_one=$1
read_two=$2
short=$3

#gunzip input/*.gz
bwa mem chr6/chr6.fa ${read_one}  ${read_two} > input/${short}.sam
echo "bwa done"
samtools view -bT chr6/chr6.fa input/${short}.sam > input/${short}.bam
samtools sort input/${short}.bam > input/${short}.sorted.bam
samtools index input/${short}.sorted.bam input/${short}.sorted.bai
echo "done index and align"

perl bin/typer.sh input/${short}.sorted.bam ${short}
