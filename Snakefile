"""Snakemake workflow documenting BIOSTAT 666 project

Note: 1K Genomes to PLINK Modified from https://www.biostars.org/p/335605/
"""

from pathlib import Path

CHROMS_INCLUDED = [str(x) for x in range(1,23)]
STUDIES = ["adhd_eur_jun2017","ocd_aug2017"]

VCF_DIR = Path("input/vcf")
ASSOC_DIR = Path("input/gwas")
INT_DIR = Path("intermediate/")
OUT_DIR = Path("output/")

MEM_LIMIT_CHROM = 1024*4
MEM_LIMIT_MERGED = 1024*64
wildcard_constraints:
        chr="\d+"
rule all:
    input: 
        expand(str( OUT_DIR / "chr{chr}_co.clumped"),chr=CHROMS_INCLUDED),
        expand(str( OUT_DIR / "chr{chr}_{study}.clumped"),chr=CHROMS_INCLUDED,study=STUDIES),

# Filter 1K genomes panel to only include european individuals
rule filter_sample_list:
    input: "input/integrated_call_samples_v3.20130502.ALL.panel"
    output: str(INT_DIR / "european_ids.txt")
    shell:
        """
        awk '{{ if ($3=="EUR") print; }}' {input} | cut -f1 > {output}
        """
rule reformat_assoc:
    input: ASSOC_DIR = ASSOC_DIR / "{study}"
    output: INT_DIR / "{study}_normed.assoc"
    shell:
        """
        awk 'NR==1; NR>1 {{$2 = $1":"$3":"$5":"$4; print}}' {input} > {output} 
        """
rule vcf_to_bcf:
    input: 
        vcf=VCF_DIR / "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf",
    output: str(INT_DIR / "raw_bcf" /  "chr{chr}.bcf")
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools view {input} -Ob > {output}
        bcftools index {output}
        """
# Rule to filter VCF to only SNVs, and only include european individuals
rule norm_filter_bcf:
    input: 
        bcf= rules.vcf_to_bcf.output,
        sample_list=rules.filter_sample_list.output
    output: INT_DIR / "normed_vcf" / "chr{chr}_normed.bcf"
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools view -S {input.sample_list} -v snps {input.bcf} | \
        bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both > {output}
        bcftools index {output}
        """

# Rule to convert BCF to PLINK format
rule vcf_to_plink:
    input: 
        vcf=rules.norm_filter_bcf.output,
    output: INT_DIR/"chr{chr}.bed"
    shell:
        """
        plink --memory {MEM_LIMIT_CHROM} --noweb \
        --bcf {input.vcf} \
        --keep-allele-order \
        --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail \
        --make-bed --out {INT_DIR}/chr{wildcards.chr}
        """

# Joint clumping of SNPs
rule co_clump_snps:
    input: 
        data=rules.vcf_to_plink.output, 
        asso=expand(str(INT_DIR / "{study}_normed.assoc"),study=STUDIES)
    output: OUT_DIR / "chr{chr}_co.clumped"
    shell:
        """
        mkdir -p $(dirname "{output}")
        plink --memory {MEM_LIMIT_CHROM} --bfile {INT_DIR}/chr{wildcards.chr} --clump {input.asso} --out {OUT_DIR}/chr{wildcards.chr}_co
        """

rule clump_snps:
    input: 
        data=rules.vcf_to_plink.output, 
        asso=INT_DIR / "{study}_normed.assoc"
    output: OUT_DIR / "chr{chr}_{study}.clumped"
    shell:
        """
        mkdir -p $(dirname "{output}")
        touch {output}
        plink --memory {MEM_LIMIT_CHROM} --bfile {INT_DIR}/chr{wildcards.chr} --clump {input.asso} --out {OUT_DIR}/chr{wildcards.chr}_{wildcards.study}
        """