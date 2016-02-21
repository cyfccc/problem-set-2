#! /usr/bin/env bash

data="../data-sets"

tfbs_bed=$data/bed/encode.tfbs.chr22.bed.gz
h3k4me3_bed=$data/bed/encode.h3k4me3.hela.chr22.bed.gz
gzcat $tfbs_bed |awk '$4 == "CTCF"'>ctcf.bed
answer_1=$(bedtools intersect -a $h3k4me3_bed -b ctcf.bed -wo\
|awk '{print $NF}'\
|sort -nr\
|head -n1) 


echo -e "chr22\t19000000\t19000500">temp_region.bed
chr22_nucleotides=$data/fasta/hg19.chr22.fa
if [[ -f $data/fasta/hg19.chr22.fa.gz ]]; then
gunzip $data/fasta/hg19.chr22.fa.gz
fi
answer_2=$(bedtools nuc -fi $chr22_nucleotides -bed temp_region.bed\
|grep -v '^#'\
|awk '{print $5}')


ctcf_signal=$data/bedtools/ctcf.hela.chr22.bg.gz
answer_3=$(bedtools map -a ctcf.bed -b $ctcf_signal -c 4 -o mean\
|sort -k5n\
|tail -n1\
|awk '{print $3-$2}')


tss=$data/bed/tss.hg19.chr22.bed.gz
genome=$data/genome/hg19.genome
gzcat $tss|awk 'BEGIN{OFS="\t"} ($6=="+"){print $1, $2, $3, $4, $6}'\
|bedtools flank -i - -g $genome -l 1000 -r 0\
|bedtools map -a - -b $ctcf_signal -c 4 -o median > temp.bed

gzcat $tss|awk 'BEGIN{OFS="\t"} ($6=="-"){print $1, $2, $3, $4, $6}'\
|bedtools flank -i - -g $genome -l 0 -r 1000\
|bedtools map -a - -b $ctcf_signal -c 4 -o median >> temp.bed
answer_4=$(cat temp.bed|sort -k6n|tail -n1|awk '{print $4}')

genes=$data/bed/genes.hg19.bed.gz
bedtools complement -i $genes -g $genome\
|awk '($1=="chr22"){print $1, $2, $3, $3-$2}'\
|sort -k4n\
|tail -n1 > temp2.bed
chr=$(cat temp2.bed|awk '{print $1}')
start=$(cat temp2.bed|awk '{print $2}')
end=$(cat temp2.bed|awk '{print $3}')


echo "answer-1: $answer_1
answer-2: $answer_2
answer-3: $answer_3
answer-4: $answer_4
answer-5: $chr:$start-$end" > answers.yml
