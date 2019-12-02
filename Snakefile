"""Snakemake workflow documenting BIOSTAT 666 project

Note: 1K Genomes to PLINK Modified from https://www.biostars.org/p/335605/
"""

from pathlib import Path

CHROMS_INCLUDED = ["21","22"]
STUDIES = ["adhd_jul2017","ocd_aug2017"]

VCF_DIR = Path("input/vcf")
ASSOC_DIR = Path("input/gwas")
INT_DIR = Path("intermediate/")
OUT_DIR = Path("output/")

# Filter 1K genomes panel to only include european individuals
rule filter_sample_list:
    input: "input/integrated_call_samples_v3.20130502.ALL.panel"
    output: str(INT_DIR / "european_ids.txt")
    shell:
        """
        awk '{{ if ($3=="EUR") print; }}' {input} | cut -f1 > {output}
        """
rule vcf_to_bcf:
    input: 
        vcf=VCF_DIR / "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
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
        bcftools view -S {input.sample_list} -v snps {input.bcf} | bcftools norm -Ob --rm-dup both > {output}
        bcftools index {output}
        """

# Rule to convert BCF to PLINK format
rule vcf_to_plink:
    input: 
        vcf=rules.norm_filter_bcf.output,
    output: INT_DIR/"chr{chr}.bed"
    shell:
        """
        plink --noweb \
        --bcf {input.vcf} --keep-allele-order \
        --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail \
        --make-bed --out {INT_DIR}/chr{wildcards.chr}
        """

# Merge PLINK chromosome files
rule merge_chroms:
    input: expand(rules.vcf_to_plink.output,chr=CHROMS_INCLUDED)
    output: INT_DIR / "merged.bed"
    shell:
        """
        rm -f {INT_DIR}/merge_list
        for line in {CHROMS_INCLUDED}
        do
            echo "{INT_DIR}/chr$line" >> {INT_DIR}/merge_list
        done
        plink --merge-list {INT_DIR}/merge_list --out {INT_DIR}/merged
        """

# Joint clumping of SNPs
rule clump_snps:
    input: data=rules.merge_chroms.output, asso=expand(str(ASSOC_DIR / "{study}"),study=STUDIES)
    output: OUT_DIR / "gwas_snps.clumped"
    shell:
        """
        mkdir -p $(dirname "{output}")
        plink --bfile {INT_DIR}/merged --clump {input.asso} --out {OUT_DIR}/gwas_snps
        """