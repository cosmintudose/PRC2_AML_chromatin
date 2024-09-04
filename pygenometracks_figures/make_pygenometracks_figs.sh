#!/bin/bash

#NREP
pyGenomeTracks --tracks ./ini_files/tracks_cut_and_run_only.ini --region chr5:111500000-112700000 -o nrep.pdf

#FLT3
pyGenomeTracks --tracks ./ini_files/tracks_cut_and_run_only.ini --region chr13:27500000-28600000 -o flt3.pdf

#SOX4
pyGenomeTracks --tracks ./ini_files/tracks_cut_and_run_only.ini --region chr6:21000000-22000000 -o sox4.pdf

#CDCA7
pyGenomeTracks --tracks ./ini_files/tracks_cut_and_run_only.ini --region chr2:172700000-173900000 -o cdca7.pdf

#SKIDA1
pyGenomeTracks --tracks ./ini_files/tracks_hic_cut_and_run_atac.ini --region chr10:21000000-22500000 -o skida1_open_frags.pdf 

#PCDH9
pyGenomeTracks --tracks ./ini_files/tracks_hic_cut_and_run_atac.ini --region chr13:61000000-68000000 -o pcdh9_open_frags.pdf  

#LIN28B
pyGenomeTracks --tracks ./ini_files/LIN28B_hic.ini --region chr6:95000000-108000000 -o lin28b.pdf