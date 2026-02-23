#create_outputs.smk

rule index_bams:
	input:
		"deduped/{sample}/{sample}_deduped.bam"
	output:
		"deduped/{sample}/{sample}_deduped.bam.bai"
	shell:
		"samtools index {input}"


rule bedgraph: #code taken from https://github.com/FredHutch/SEACR
	input:
		"deduped/{sample}/{sample}_deduped.bam"
	output:
		"bedgraphs/{sample}/{sample}.bedgraph"
	params:
		genome = "/home/administrator/Bowtie2Index/chrom.sizes"
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
		genome = "/home/administrator/Bowtie2Index/chrom.sizes"
	shell:
		"bedGraphToBigWig {input} "
		"{params.genome} "
		"{output} "

rule bamcov:
	input:
		input_files = "deduped/{sample}/{sample}_deduped.bam",
	output:
		"bigwigs/{sample}.bigwig"
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
		input_files = "bedgraphs/{condition}_{mark}/{condition}_{mark}.bedgraph",

		input_control = "bedgraphs/{condition}_IgG_1/{condition}_IgG_1.bedgraph"
	output:
		"seacr/{condition}_{mark}"
	params:
		prefix = ".stringent.bed"
	shell:
		"/home/administrator/SEACR-master/SEACR_1.3.sh {input.input_files} "
		"{input.input_control} "
		"norm stringent "
		"{output}{params.prefix}"
