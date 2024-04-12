configfile: "config.yaml"

# rule all:
#     input:
#         expand("{output_dir}/{sample}.fa.gz", output_dir=config['output_dir'], sample=config['samples'])

rule all:
    input:
        expand("results/denovo.haplotagged.{sample}_mhc.500k.bam", sample=config['samples'])

# Filter reads from .fire.bam file to only those within selected region.
# samtools view -@ `nproc` -b ../data/fire-bams/GM12878.fire.bam -ML ../data/mhc.500k.slop.bed > example.bam
rule extract_bam:
    input:
        bam=lambda wildcards: config['input_file_dir'] + "/" + wildcards.sample + ".fire.bam",
        bed=config['genome_region']
    output:
        bam=temp("{working_dir}/{sample}.bam")
    conda:
        config['conda_env']
    shell:
        """
        samtools view -@ $(nproc) -b {input.bam} -L {input.bed} > {output.bam}
        """

# Generate a compressed fasta file from the bam file.
# samtools fasta example.bam | bgzip -@ `nproc` > example.fa.gz
rule convert_to_fasta:
    input:
        bam="{working_dir}/{sample}.bam"
    output:
        fasta_gz="{working_dir}/{sample}.fa.gz"
    conda:
        config['conda_env']
    shell:
        """
        samtools fasta {input.bam} | bgzip -@ $(nproc) > {output.fasta_gz}
        """

# Run hifiasm on filtered reads.
# mkdir -p results
# hifiasm -o results/example -t `nproc` example.fa.gz
rule run_hifiasm:
    input:
        fasta="{working_dir}/{sample}.fa"
    output:
        asm="{output_dir}/{sample}.asm"
    conda:
        config['conda_env']
    shell:
        """
        mkdir -p {output_dir}
        hifiasm -o {output.asm} -t $(nproc) {input.fasta}
        """

# rule compress_fasta:
#     input:
#         fasta="{working_dir}/{sample}.fa"
#     output:
#         fasta_gz="{output_dir}/{sample}.fa.gz"
#     conda:
#         "env.yml"
#     shell:
#         """
#         bgzip -@ $(nproc) {input.fasta} > {output.fasta_gz}
#         """

# Index .fire.bam file, if needed.
rule index_bam:
    input:
        "{input_file_dir}/{sample}.fire.bam"
    output:
        "{input_file_dir}/{sample}.fire.bam.bai"
    conda:
        config['conda_env']
    shell:
        "samtools index {input}"


# map reads back to the denovo assembly
# pbmm2 align --unmapped --sort -j `nproc` denovo-mhc.fasta example.bam denovo.example.bam 
rule map_reads_to_assembly:
    input:
        fasta="results/denovo-{sample}_mhc.500k.fasta",
        bam="{input_file_dir}/{sample}_mhc.500k.bam"
    output:
        aligned_bam="results/denovo.{sample}_mhc.500k.bam"
    conda:
        config['conda_env']
    shell:
        """
        pbmm2 align --unmapped --sort -j $(nproc) {input.fasta} {input.bam} {output.aligned_bam}
        """

# Split reads that align with our region by haplotype, and add them to a fasta file.
# gfatools gfa2fa results/example.bp.hap1.p_ctg.gfa  > denovo-mhc.fasta
# gfatools gfa2fa results/example.bp.hap2.p_ctg.gfa >> denovo-mhc.fasta
rule gfa_to_fasta:
    input:
        hap1="results/{sample}_mhc.500k.bp.hap1.p_ctg.gfa",
        hap2="results/{sample}_mhc.500k.bp.hap2.p_ctg.gfa"
    output:
        fasta="results/denovo-{sample}_mhc.500k.fasta"
    conda:
        config['conda_env']
    shell:
        """
        gfatools gfa2fa {input.hap1} > {output.fasta}
        gfatools gfa2fa {input.hap2} >> {output.fasta}
        """

rule tag_reads_by_haplotype:
    input:
        bam="results/denovo.{sample}_mhc.500k.bam"
    output:
        tagged_bam="results/denovo.haplotagged.{sample}_mhc.500k.bam"
    conda:
        config['conda_env']
    shell:
        """
        scripts/haplotag-reads-by-asm.py --hap1-tag h1tg --hap2-tag h2tg {input.bam} {output.tagged_bam}
        """


