#!/bin/bash
 
export PATH=/home/fgao/software/STAR-2.7.8a/bin/Linux_x86_64/:$PATH

# star index
mkdir -p /home/fgao/genome_reference/gencode_grcm39/STARIndex_50
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /home/fgao/genome_reference/gencode_grcm39/STARIndex_50 \
--genomeFastaFiles /home/fgao/genome_reference/gencode_grcm39/GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile /home/fgao/genome_reference/gencode_grcm39/gencode.vM28.annotation.gtf \
--sjdbOverhang 50

mkdir -p fastq_trim
# star alignment
while read line;
do
ext1="_R1_001.fastq.gz"
file1=$line$ext1
zcat $file1 | awk '{if(NR%2==0) print substr($1,1,50); else print}' | gzip > fastq_trim/$file1
done < sample_ID.txt

while read line;
do
ext1="_R1_001.fastq.gz"
file1=$line$ext1
STAR --runThreadN 16 \
--genomeDir /home/fgao/genome_reference/gencode_grcm39/STARIndex_50 \
--readFilesIn fastq_trim/$file1 \
--readFilesCommand zcat \
--outFileNamePrefix $line \
--outSAMtype BAM SortedByCoordinate

#2nd pass
STAR --runThreadN 16 \
--genomeDir /home/fgao/genome_reference/gencode_grcm39/STARIndex_50 \
--quantMode TranscriptomeSAM \
--sjdbFileChrStartEnd ${line}SJ.out.tab \
--readFilesCommand zcat \
--outFileNamePrefix ${line}_p2 \
--readFilesIn fastq_trim/$file1 \
--outSAMtype BAM SortedByCoordinate

samtools index ${line}Aligned.sortedByCoord.out.bam
samtools index ${line}_p2Aligned.sortedByCoord.out.bam
done < sample_ID.txt

#rsem indexing
awk '{if($3=="transcript") print $10"\t"$12}' /home/fgao/genome_reference/gencode_grcm39/gencode.vM28.annotation.gtf | sed 's/"//g' - | sed 's/;//g' - > knownIsoforms.txt

/home/fgao/software/RSEM-1.3.3/rsem-prepare-reference --gtf /home/fgao/genome_reference/gencode_grcm39/gencode.vM28.annotation.gtf \
                            --transcript-to-gene-map knownIsoforms.txt \
                            --star \
                            --star-path /home/fgao/software/STAR-2.7.8a/bin/Linux_x86_64/ \
                            -p 8 \
                            /home/fgao/genome_reference/gencode_grcm39/GRCm39.primary_assembly.genome.fa \
                            /home/fgao/genome_reference/gencode_grcm39/rsem_index

#rsem gene quantification
while read line;
do
/home/fgao/software/RSEM-1.3.3/rsem-calculate-expression -p 32 \
        --alignments ${line}_p2Aligned.toTranscriptome.out.bam /home/fgao/genome_reference/gencode_grcm39/rsem_index $line
done < sample_ID.txt

awk '{if($3=="gene") print $10"\t"$14}' /home/fgao/genome_reference/gencode_grcm39/gencode.vM28.annotation.gtf | sed 's/"//g' - | sed 's/;//g' - > gene_id_name.txt
