# PROVENANCE: these 4 comment lines were added in this snapshot and are not part of the original.
# KaryoScope-archive workflow/rules/pangenome.smk at commit 0bf3a7d ("master pangenome QC table").
# Original sha256 286ca01346c0301d417859f47ccc3b3aff7f67ea1a3c5eaa938fc14e91d16df0; drop the first 4 lines to get it.
# Kept as a record; build_master_qc_table.py is the runnable version.
HAPS = ["hap1", "hap2"]
TOOLS = ["flagger_hifi", "flagger_ont", "nucflag"]
REGIONS = ["subtelo_p", "centromere", "subtelo_q"]

rule get_error_free_T2T_contigs:
  input:
    pangenome_fais,
    get_all_sample_QCs
  output:
    "results/pangenome/error_free/{qc}_error_free_T2T_contigs.txt"
  params:
    samples=pangenome_samples
  script:
    "../scripts/get_error_free_T2T_contigs.py"

rule get_regions:
  input:
    "results/{sample}/KaryoScope/human_assembly_haploid/{sample}.human_assembly_haploid.binsize1000000.scaffold_stats.tsv",
    "results/{sample}/assembly/{hap}/{sample}.{hap}.fasta.gz.fai",
    "results/{sample}/KaryoScope/human_assembly_haploid/{sample}.human_assembly_haploid.featureIDs.binsize1000000.centromeres.scaffolded.bed"
  output:
    "results/pangenome/errors/{sample}/{hap}/{chrom}/{sample}.{hap}.{chrom}.subtelo_p.bed",
    "results/pangenome/errors/{sample}/{hap}/{chrom}/{sample}.{hap}.{chrom}.subtelo_q.bed",
    "results/pangenome/errors/{sample}/{hap}/{chrom}/{sample}.{hap}.{chrom}.centromere.bed"
  resources:
    runtime=5
  script:
    "../scripts/get_regions.py"

rule find_errors:
  input:
    a=get_QC_wildcards,
    b="results/pangenome/errors/{sample}/{hap}/{chrom}/{sample}.{hap}.{chrom}.{zoom}.bed"
  output:
    "results/pangenome/errors/{sample}/{hap}/{chrom}/{sample}.{hap}.{chrom}.{zoom}.{qc}.errors.bed",
  resources:
    runtime=15
  shell:
    "bedtools intersect -a {input.a} -b {input.b} > {output}"

rule interesting_example_errors:
  input:
    interesting_example_errors
  output:
    "results/pangenome/errors/interesting/interesting.centromere.nucflag.errors.bed"
  script:
    "../scripts/interesting_example_errors.py"

rule interesting_example_errors_subtelomere:
  input:
    interesting_example_errors_subtelomere
  output:
    "results/pangenome/errors/interesting/interesting.subtelomere.nucflag.errors.bed"
  script:
    "../scripts/interesting_example_errors.py"

rule aggregate_master_qc_table:
  input:
    # 1. All intersection files from the find_errors rule
    intersections=expand(
        "results/pangenome/errors/{sample}/{hap}/{chrom}/{sample}.{hap}.{chrom}.{zoom}.{qc}.errors.bed",
        sample=pangenome_samples,
        hap=HAPS,
        chrom=CHROMOSOME_LIST,
        zoom=REGIONS,
        qc=TOOLS
    ),
    # 2. All scaffold stats files to build the T2T map
    scaffold_stats=expand(
        "results/{sample}/KaryoScope/human_assembly_haploid/{sample}.human_assembly_haploid.binsize1000000.scaffold_stats.tsv",
        sample=pangenome_samples
    )
  output:
    "results/pangenome/master_qc_table.full.tsv"
  params:
    samples=pangenome_samples,
    haps=HAPS,
    chroms=CHROMOSOME_LIST,
    tools=TOOLS,
    regions=REGIONS
  script:
    "../scripts/aggregate_master_qc.py"
