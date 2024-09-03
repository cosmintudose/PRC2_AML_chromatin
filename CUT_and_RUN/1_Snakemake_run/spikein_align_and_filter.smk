#spikein_align_and_filter.smk

rule align_reads_spikein:
	input:
		read1 = "fastq_trimmed/{sample}_R1_paired.fq.gz",
		read2 = "fastq_trimmed/{sample}_R2_paired.fq.gz",
		wait_for_finish = rules.align_reads.output
	output:
		alignment = temp("sams_spikein/spikein_{sample}.sam"), #removes temp sam after the .bams are created
		stats = "stats_spikein_{sample}_bowtie2.txt"
	threads: 8
	params:
		yeast_genome_idx = "/mnt/data/yeast_genome_bowtie2/SacCer" #path to genome sample
	shell:
		"bowtie2 --end-to-end --dovetail -I 10 -X 700 --no-mixed --no-discordant --threads {threads} -x {params.yeast_genome_idx} -1 {input.read1} -2 {input.read2} "
		"-S {output.alignment} &> {output.stats}"


rule sam_to_bam_spikein:
	input:
		rules.align_reads_spikein.output.alignment
	output:
		"bams_spikein/spikein_{sample}.bam"
	shell:
		"samtools view -bS {input} > {output}"

rule count_reads_after_alignment_spikein: 
	input:
		"bams_spikein/spikein_{sample}.bam"
	output:
		"count_reads/after_alignment/spikein/{sample}_reads_after_alignment.txt"
	shell:
		"bamtools stats -in {input} > {output}"



rule sort_bam_spikein:
    input:
        bam_files = "bams_spikein/spikein_{sample}.bam"
    output:
        temp("sorted_bams_spikein/{sample}/{sample}_sorted.bam")
    shell:
        "java -jar /mnt/data/scripts/picard.jar SortSam I={input.bam_files} "
        "O={output} "
        "SORT_ORDER=coordinate "
        "VALIDATION_STRINGENCY=LENIENT"


rule mark_duplicates_spikein:
	input:
		rules.sort_bam_spikein.output
	output:
		dupes = "dupes_spikein/{sample}/{sample}_marked_duplicates.bam",
		metrics = "dupes_spikein/{sample}_dup_metrics.txt"
	shell:
		"java -jar /mnt/data/scripts/picard.jar MarkDuplicates I={input} "
		"O={output.dupes} "
		"M={output.metrics} "
		"VALIDATION_STRINGENCY=LENIENT"


rule index_bams_dupes_spikein:
	input:
		rules.mark_duplicates_spikein.output.dupes
	output:
		"dupes_spikein/{sample}/{sample}_marked_duplicates.bam.bai"
	shell:
		"samtools index {input}"


rule filter_duplicates_spikein:
	input:
		input_files = "dupes_spikein/{sample}/{sample}_marked_duplicates.bam",
		index_files = "dupes_spikein/{sample}/{sample}_marked_duplicates.bam.bai"
	output:
		"deduped_spikein/{sample}/{sample}_deduped.bam"
	shell:
		"samtools view -F 1280 -f 2 -b {input.input_files} > {output}" #remove multimappers -F 1280 (multimappers & pcr/optical duplicates) -F 3328 remove supplementary alignment as well


rule index_bams_spikein:
	input:
		rules.filter_duplicates_spikein.output
	output:
		"deduped_spikein/{sample}/{sample}_deduped.bam.bai"
	shell:
		"samtools index {input}"


rule count_reads_after_removing_dupes_spikein: 
	input:
		input_files = rules.filter_duplicates_spikein.output,
		wait_for_finish = rules.index_bams_spikein.output
	output:
		"count_reads/after_alignment/spikein/{sample}_reads_deduped_bam.txt"
	shell:
		"bamtools stats -in {input.input_files} > {output}"