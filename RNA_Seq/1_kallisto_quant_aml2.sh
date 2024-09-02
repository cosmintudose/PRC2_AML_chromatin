cd ./RNA_Seq

#checks quality of fastq files
fastqc *.gz -t 4 #number of threads it's using, can be changed
echo "fastqc finished"

#builds index for reference genome
kallisto index -i Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa #change name if a different version, do the same downstream
echo "Index built"

#pseudoalignment and output in separate directories - this version also creates .bam files, but they're not necessary
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./A1_A2 -t 4 ./A1_A2_1.fq.gz ./A1_A2_2.fq.gz &> ./A1_A2.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./A3_A4 -t 4 ./A3_A4_1.fq.gz ./A3_A4_2.fq.gz &> ./A3_A4.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./A5_A6 -t 4 ./A5_A6_1.fq.gz ./A5_A6_2.fq.gz &> ./A5_A6.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./B1_B2 -t 4 ./B1_B2_1.fq.gz ./B1_B2_2.fq.gz &> ./B1_B2.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./B3_B4 -t 4 ./B3_B4_1.fq.gz ./B3_B4_2.fq.gz &> ./B3_B4.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./B5_B6 -t 4 ./B5_B6_1.fq.gz ./B5_B6_2.fq.gz &> ./B5_B6.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./C1_C2 -t 4 ./C1_C2_1.fq.gz ./C1_C2_2.fq.gz &> ./C1_C2.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./C3_C4 -t 4 ./C3_C4_1.fq.gz ./C3_C4_2.fq.gz &> ./C3_C4.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o ./C5_C6 -t 4 ./C5_C6_1.fq.gz ./C5_C6_2.fq.gz &> ./C5_C6.log

echo "Alignment finished"

#summarizes qc and mapping results in single html
multiqc -d .

