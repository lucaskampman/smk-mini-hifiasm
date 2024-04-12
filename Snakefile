configfile: "config.yaml"

# rule all:
#     input:
#         expand("{output_dir}/{sample}.fa.gz", output_dir=config['output_dir'], sample=config['samples'])

rule all:
    input:
        expand("results/denovo.haplotagged.{sample}_mhc.500k.bam", sample=config['samples'])

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

rule index_bam:
    input:
        "{input_file_dir}/{sample}.fire.bam"
    output:
        "{input_file_dir}/{sample}.fire.bam.bai"
    conda:
        config['conda_env']
    shell:
        "samtools index {input}"

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
        haplotag-reads-by-asm.py --hap1-tag h1tg --hap2-tag h2tg {input.bam} {output.tagged_bam}
        """


