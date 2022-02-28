#!/bin/bash
 
export PATH=/home/fgao/software/STAR-2.7.8a/bin/Linux_x86_64/:$PATH

# star index
mkdir -p /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/STARIndex_50
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/STARIndex_50 \
--genomeFastaFiles /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
--sjdbGTFfile /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf \
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
--genomeDir /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/STARIndex_50 \
--readFilesIn fastq_trim/$file1 \
--readFilesCommand zcat \
--outFileNamePrefix ${line}_mm10 \
--outSAMtype BAM SortedByCoordinate

#2nd pass
STAR --runThreadN 16 \
--genomeDir /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/STARIndex_50 \
--quantMode TranscriptomeSAM \
--sjdbFileChrStartEnd ${line}_mm10SJ.out.tab \
--readFilesCommand zcat \
--outFileNamePrefix ${line}_mm10p2 \
--readFilesIn fastq_trim/$file1 \
--outSAMtype BAM SortedByCoordinate

samtools index ${line}_mm10Aligned.sortedByCoord.out.bam
samtools index ${line}_mm10p2Aligned.sortedByCoord.out.bam
done < sample_ID.txt

#rsem indexing
awk '{if($0~"gene_id" && $0~"transcript_id") print $10"\t"$(NF-2)}' /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf | sort | uniq | sed
 's/"//g' - | sed 's/;//g' - > knownIsoforms_mm10.txt
/home/fgao/software/RSEM-1.3.3/rsem-prepare-reference --gtf /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf \
                            --transcript-to-gene-map knownIsoforms_mm10.txt \
                            --star \
                            --star-path /home/fgao/software/STAR-2.7.8a/bin/Linux_x86_64/ \
                            -p 8 \
                            /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
                            /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/rsem_index

#rsem gene quantification
while read line;
do
/home/fgao/software/RSEM-1.3.3/rsem-calculate-expression -p 32 \
        --alignments ${line}_mm10p2Aligned.toTranscriptome.out.bam /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/rsem_index ${line}_mm10
done < sample_ID.txt
