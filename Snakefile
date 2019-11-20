"""Snakemake workflow documenting BIOSTAT 666 project
"""

from pathlib import Path

CHROMS_INCLUDED = ["21","22"]
STUDIES = ["adhd_jul2017","ocd_aug2017"]

VCF_DIR = Path("input/vcf")
ASSOC_DIR = Path("input/gwas")
INT_DIR = Path("intermediate/")
REF_FASTA = Path("input/human_g1k_v37.fasta")
OUT_DIR = Path("output/")
rule norm_vcf:
    input: vcf=VCF_DIR / "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",ref=REF_FASTA
    output: INT_DIR / "chr{chr}_normed.bcf"
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools norm -m-any --check-ref w -f {input.ref} {input.vcf} | \
        bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both > {output}
        bcftools index {output}
        """
rule vcf_to_plink:
    input: rules.norm_vcf.output
    output: INT_DIR/"chr{chr}.bed"
    shell:
        """
        plink --noweb --bcf {input} \
         --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed \
        --out {INT_DIR}/chr{wildcards.chr}
        """
rule reformat_assoc:
    input: ASSOC_DIR = ASSOC_DIR / "{study}"
    output: INT_DIR / "{study}_normed.assoc"
    shell:
        """
        awk 'NR==1; NR>1 {{$2 = $1":"$3":"$5":"$4; print}}' {input} > {output} 
        """
rule merge_chroms:
    input: expand(rules.vcf_to_plink.output,chr=CHROMS_INCLUDED)
    output: INT_DIR / "merged.bed"
    shell:
        """
        for line in {CHROMS_INCLUDED}
        do
            echo "{INT_DIR}/chr$line \n" >> {INT_DIR}/merge_list
        done
        plink --merge-list {INT_DIR}/merge_list --out {INT_DIR}/merged
        """
rule clump_snps:
    input: data=rules.merge_chroms.output, asso=expand(rules.reformat_assoc.output,study=STUDIES)
    output: OUT_DIR / "gwas_snps.clumped"
    shell:
        """
        mkdir -p $(dirname "{output}")
        plink --bfile {INT_DIR}/merged --clump {input.asso} --out {OUT_DIR}/gwas_snps
        """