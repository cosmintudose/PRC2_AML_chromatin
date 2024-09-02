#creating directories needed

mkdir ./RNA_Seq/results_files/
mkdir ./RNA_Seq/plots/
mkdir publicly_available_data

#MOLM13 CUT&RUN
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6893nnn/GSM6893214/suppl/GSM6893214%5FH3K27Tri%2D9733S.SeqDepthNorm.bw --directory-prefix=./publicly_available_data
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6893nnn/GSM6893214/suppl/GSM6893214%5FH3K27Tri%2D9733S%5FInput%2Dq05%5Fpeaks.broadPeak.gz --directory-prefix=./publicly_available_data

#HL60 chip-seq
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5330nnn/GSM5330834/suppl/GSM5330834%5FENCFF598BHR%5Fpseudoreplicated%5Fpeaks%5FGRCh38.bigBed --directory-prefix=./publicly_available_data

#K562 Micro-C - this is 7.6 GB
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206131/suppl/GSE206131%5FK562.mapq%5F30.mcool --directory-prefix=./publicly_available_data

#TARGET AML data
wget https://cbioportal-datahub.s3.amazonaws.com/aml_target_2018_pub.tar.gz --directory-prefix=./publicly_available_data

tar -xvzf ./publicly_available_data/aml_target_2018_pub.tar.gz -C ./publicly_available_data