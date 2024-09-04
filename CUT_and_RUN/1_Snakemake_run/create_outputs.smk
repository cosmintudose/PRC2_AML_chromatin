#create_outputs.smk

rule bedgraph: #code taken from https://github.com/FredHutch/SEACR
	input:
		"deduped/{sample}/{sample}_deduped.bam"
	output:
		"bedgraphs/{sample}/{sample}.bedgraph"
	params:
		genome = "../../genome/chrom.sizes"
	shell:
		"samtools sort -n {input} | "
		"bedtools bamtobed -bedpe -i | "
        "awk '$1==$4 && $6-$2 < 1000 {{print $0}}' | "
        "cut -f 1,2,6 | "
        "sort -k1,1 -k2,2n -k3,3n | "
        "bedtools genomecov -bg -i - -g {params.genome} > {output}"

rule bedgraph_to_bigwig: 
	input: 
		rules.bedgraph.output
	output:
		"bedgraphs/{sample}/{sample}.bigwig"
	params:
		genome = "../../genome/chrom.sizes"
	shell:
		"bedGraphToBigWig {input} "
		"{params.genome} "
		"{output} "

rule bamcov:
	input:
		input_files = "deduped/{sample}/{sample}_deduped.bam",
		wait_for_finish = "count_reads/after_alignment/{sample}_reads_deduped_bam.txt"
	output:
		"bigwigs/{sample}_deduped.bigwig"
	shell:
		"bamCoverage --bam {input.input_files} "
		"-o {output} " 
		"--binSize 30 "
		"--smoothLength 60 "
		"--normalizeUsing CPM "
		"--effectiveGenomeSize 2913022398 "
		"--extendReads"


rule seacr:
	input:
		input_files = "bedgraphs/{condition}_{mark}_1/{condition}_{mark}_1.bedgraph",
		input_control = "bedgraphs/{condition}_igg_ctrl_1/{condition}_igg_ctrl_1.bedgraph"
	output:
		"seacr/{condition}_{mark}"
	shell:
		"./SEACR_1.3.sh {input.input_files} "
		"{input.input_control} "
		"norm stringent "
		"{output}"
