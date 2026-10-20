#repeatmasker
rule get_repeatmasker_HG002_hap1:
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.hap1.RepeatMasker.bb"
  resources:
    mem_mb=5000,
    runtime=5
  shell:
    "wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/repeatmasker/hg002v1.1.paternal.bb -O {output}"

rule get_repeatmasker_HG002_hap2:
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.hap2.RepeatMasker.bb"
  resources:
    mem_mb=5000,
    runtime=5
  shell:
    "wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/repeatmasker/hg002v1.1.maternal.bb -O {output}"

rule repeatmasker_HG002_bb_to_bed:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.bb"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.bed"
  resources:
    mem_mb=5000,
    runtime=10
  shell:
    "bigBedToBed {input} {output}"

rule process_repeatmasker_HG002:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.bed"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.bed"
  shell:
    "python {workflow.basedir}/workflow/scripts/process_repeatmasker.py --input {input} --output {output}"

rule sort_repeatmasker_HG002:
  input:
    bed=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.bed",
    fai=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.fasta.gz.fai"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.sorted.bed"
  shell:
    "bedtools sort -i {input.bed} -g {input.fai} > {output}"

rule merge_repeatmasker_HG002:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.sorted.bed"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.sorted.merged.bed"
  shell:
    "python {workflow.basedir}/workflow/scripts/priority_merge.py --input {input} --order rRNA,scRNA,snRNA,srpRNA,tRNA,Satellite,RC,Retroposon,DNA,LTR,SINE,LINE,Simple_repeat,Low_complexity,Unknown --output {output}"

rule find_n_regions_HG002:
  input:
    "results/HG002/assembly/{hap}/HG002.{hap}.fasta.gz"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.n_regions.bed"
  shell:
    "python {workflow.basedir}/workflow/scripts/find_n_regions.py {input} > {output}"

rule fill_gaps_HG002:
  input:
    fai=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.fasta.gz.fai",
    bed=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.sorted.merged.bed",
    gaps=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.n_regions.bed"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.sorted.merged.filled.bed"
  params:
    kmer_size=config.get("kmer_size", 31)
  shell:
    "python {workflow.basedir}/workflow/scripts/fill_gaps.py --fai {input.fai} --input {input.bed} --output {output} --kmer-size {params.kmer_size} --gaps {input.gaps}"

rule compare_HG002:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.RepeatMasker.processed.sorted.merged.filled.bed",
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.repeat.presmoothed.features.bed.gz",
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.repeat.smoothed.features.bed.gz"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.bed"
  shell:
    "bedtools unionbedg -i {input} > {output}"

rule contrast_HG002:
  input:
    bed="results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.bed",
    hierarchy=Path(workflow.basedir)/"resources/databases/KS_human_CHM13_v2/KS_human_CHM13_v2.hierarchy.txt"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.contrast.tsv"
  params:
    merge="Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat"
  shell:
    "python {workflow.basedir}/workflow/scripts/contrast_HG002.py --input {input.bed} --hierarchy {input.hierarchy} --feature-set repeat --merge '{params.merge}' --output {output}"

#CenSat

rule get_censat_HG002:
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.censat.bb"
  resources:
    mem_mb=5000,
    runtime=5
  shell:
    "wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.noheader.bb -O {output}"

rule censat_HG002_bb_to_bed:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.censat.bb"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.censat.bed"
  resources:
    mem_mb=5000,
    runtime=10
  shell:
    "bigBedToBed {input} {output}"

rule get_censat_HG002_hap1:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.censat.bed"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.hap1.censat.bed"
  shell:
    """
    grep "_PATERNAL" {input} > {output}
    """

rule get_censat_HG002_hap2:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.censat.bed"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.hap2.censat.bed"
  shell:
    """
    grep "_MATERNAL" {input} > {output}
    """

rule sort_censat_HG002:
  input:
    bed=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.bed",
    fai=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.fasta.gz.fai"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.sorted.bed"
  resources:
    mem_mb=5000,
    runtime=5
  shell:
    "bedtools sort -i {input.bed} -g {input.fai} > {output}"

rule process_censat_HG002:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.sorted.bed"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.processed.sorted.bed"
  resources:
    mem_mb=5000,
    runtime=5
  shell:
    "python {workflow.basedir}/workflow/scripts/process_censat.py --input {input} --output {output}"

rule fill_gaps_censat_HG002:
  input:
    fai=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.fasta.gz.fai",
    bed=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.processed.sorted.bed",
    gaps=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.n_regions.bed"
  output:
    filled=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.processed.sorted.filled.bed",
    centromeres=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.centromeres.bed"
  params:
    kmer_size=config.get("kmer_size", 31)
  shell:
    "python {workflow.basedir}/workflow/scripts/fill_gaps_censat.py --fai {input.fai} --input {input.bed} --output {output.filled} --centromeres {output.centromeres} --kmer-size {params.kmer_size} --gaps {input.gaps}"

rule censat_centromeres:
  input:
    censat=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.processed.sorted.filled.bed",
    centromeres=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.centromeres.bed"
  output:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.processed.sorted.filled.centromeres.bed"
  shell:
    "bedtools intersect -a {input.censat} -b {input.centromeres} > {output}"

rule presmoothed_centromeres:
  input:
    presmoothed="results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.presmoothed.features.bed.gz",
    centromeres=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.centromeres.bed"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.presmoothed.features.centromeres.bed"
  shell:
    "bedtools intersect -a {input.presmoothed} -b {input.centromeres} > {output}"

rule smoothed_centromeres:
  input:
    smoothed="results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.smoothed.features.bed.gz",
    centromeres=Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.centromeres.bed"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.smoothed.features.centromeres.bed"
  shell:
    "bedtools intersect -a {input.smoothed} -b {input.centromeres} > {output}"

rule compare_censat_HG002:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.processed.sorted.filled.bed",
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.presmoothed.features.bed.gz",
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.smoothed.features.bed.gz"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.bed"
  shell:
    "bedtools unionbedg -i {input} > {output}"

rule contrast_censat_HG002:
  input:
    bed="results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.bed",
    hierarchy=Path(workflow.basedir)/"resources/databases/KS_human_CHM13_v2/KS_human_CHM13_v2.hierarchy.txt"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.contrast.tsv"
  shell:
    "python {workflow.basedir}/workflow/scripts/contrast_HG002.py --input {input.bed} --hierarchy {input.hierarchy} --feature-set region --output {output}"

rule compare_censat_HG002_centromeres:
  input:
    Path(workflow.basedir)/"resources/KS_human_HG002/HG002.{hap}.censat.processed.sorted.filled.centromeres.bed",
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.presmoothed.features.centromeres.bed",
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.region.smoothed.features.centromeres.bed"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.centromeres.bed"
  shell:
    "bedtools unionbedg -i {input} > {output}"

rule contrast_censat_HG002_centromeres:
  input:
    bed="results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.centromeres.bed",
    hierarchy=Path(workflow.basedir)/"resources/databases/KS_human_CHM13_v2/KS_human_CHM13_v2.hierarchy.txt"
  output:
    "results/HG002/assembly/{hap}/KaryoScope/KS_human_CHM13_v2/HG002.{hap}.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.contrast.centromeres.tsv"
  shell:
    "python {workflow.basedir}/workflow/scripts/contrast_HG002.py --input {input.bed} --hierarchy {input.hierarchy} --feature-set region --output {output}"

#Figures

rule figure2_benchmarking:
  input:
    rm_hap1="results/HG002/assembly/hap1/KaryoScope/KS_human_CHM13_v2/HG002.hap1.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.contrast.tsv",
    rm_hap2="results/HG002/assembly/hap2/KaryoScope/KS_human_CHM13_v2/HG002.hap2.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.contrast.tsv",
    cs_hap1="results/HG002/assembly/hap1/KaryoScope/KS_human_CHM13_v2/HG002.hap1.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.contrast.tsv",
    cs_hap2="results/HG002/assembly/hap2/KaryoScope/KS_human_CHM13_v2/HG002.hap2.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.contrast.tsv",
    hierarchy=Path(workflow.basedir)/"resources/databases/KS_human_CHM13_v2/KS_human_CHM13_v2.hierarchy.txt"
  output:
    pdf="results/HG002/KaryoScope/KS_human_CHM13_v2/figures/figure2_benchmarking.pdf",
    png="results/HG002/KaryoScope/KS_human_CHM13_v2/figures/figure2_benchmarking.png",
    panel_a="results/HG002/KaryoScope/KS_human_CHM13_v2/figures/figure2A_repeatmasker_concordance.pdf",
    panel_b="results/HG002/KaryoScope/KS_human_CHM13_v2/figures/figure2B_censat_concordance.pdf",
    panel_c="results/HG002/KaryoScope/KS_human_CHM13_v2/figures/figure2C_confusion_repeatmasker.pdf",
    panel_d="results/HG002/KaryoScope/KS_human_CHM13_v2/figures/figure2D_confusion_censat.pdf",
    supp="results/HG002/KaryoScope/KS_human_CHM13_v2/figures/supplementary_smoothing_comparison.pdf"
  params:
    outdir="results/HG002/KaryoScope/KS_human_CHM13_v2/figures"
  shell:
    """
    Rscript {workflow.basedir}/workflow/scripts/figure2_benchmarking.R \
      --rm_hap1 {input.rm_hap1} \
      --rm_hap2 {input.rm_hap2} \
      --cs_hap1 {input.cs_hap1} \
      --cs_hap2 {input.cs_hap2} \
      --hierarchy {input.hierarchy} \
      --outdir {params.outdir}
    """

# ============================================================================
# CHM13 self-benchmark (sanity check: database vs. its own annotations)
# ============================================================================

rule trim_repeat_CHM13:
  input:
    bed=Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13_rep.bed",
    fai=Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13.fasta.gz.fai"
  output:
    Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13_rep.trimmed.bed"
  params:
    kmer_size=config.get("kmer_size", 31)
  shell:
    "python {workflow.basedir}/workflow/scripts/trim_bed_ends.py --input {input.bed} --fai {input.fai} --kmer-size {params.kmer_size} --output {output}"

rule trim_region_CHM13:
  input:
    bed=Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13_reg.bed",
    fai=Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13.fasta.gz.fai"
  output:
    Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13_reg.trimmed.bed"
  params:
    kmer_size=config.get("kmer_size", 31)
  shell:
    "python {workflow.basedir}/workflow/scripts/trim_bed_ends.py --input {input.bed} --fai {input.fai} --kmer-size {params.kmer_size} --output {output}"

rule compare_CHM13:
  input:
    Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13_rep.trimmed.bed",
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.repeat.presmoothed.features.bed.gz",
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.repeat.smoothed.features.bed.gz"
  output:
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.bed"
  shell:
    "bedtools unionbedg -i {input} > {output}"

rule contrast_CHM13:
  input:
    bed="results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.bed",
    hierarchy=Path(workflow.basedir)/"resources/databases/KS_human_CHM13_v2/KS_human_CHM13_v2.hierarchy.txt"
  output:
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.contrast.tsv"
  params:
    merge="Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat"
  shell:
    "python {workflow.basedir}/workflow/scripts/contrast_HG002.py --input {input.bed} --hierarchy {input.hierarchy} --feature-set repeat --merge '{params.merge}' --output {output}"

rule compare_censat_CHM13:
  input:
    Path(workflow.basedir)/"resources/KS_human_CHM13/CHM13_reg.trimmed.bed",
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.region.presmoothed.features.bed.gz",
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.region.smoothed.features.bed.gz"
  output:
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.bed"
  shell:
    "bedtools unionbedg -i {input} > {output}"

rule contrast_censat_CHM13:
  input:
    bed="results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.bed",
    hierarchy=Path(workflow.basedir)/"resources/databases/KS_human_CHM13_v2/KS_human_CHM13_v2.hierarchy.txt"
  output:
    "results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.contrast.tsv"
  shell:
    "python {workflow.basedir}/workflow/scripts/contrast_HG002.py --input {input.bed} --hierarchy {input.hierarchy} --feature-set region --output {output}"

rule supplementary_CHM13_benchmarking:
  input:
    rm="results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.repeatmasker.vs.presmoothed.vs.smoothed.contrast.tsv",
    cs="results/CHM13/assembly/hap1/KaryoScope/KS_human_CHM13_v2/CHM13.hap1.KS_human_CHM13_v2.censat.vs.presmoothed.vs.smoothed.contrast.tsv",
    hierarchy=Path(workflow.basedir)/"resources/databases/KS_human_CHM13_v2/KS_human_CHM13_v2.hierarchy.txt"
  output:
    pdf="results/CHM13/KaryoScope/KS_human_CHM13_v2/figures/supplementary_CHM13_benchmarking.pdf",
    png="results/CHM13/KaryoScope/KS_human_CHM13_v2/figures/supplementary_CHM13_benchmarking.png"
  params:
    outdir="results/CHM13/KaryoScope/KS_human_CHM13_v2/figures"
  shell:
    """
    Rscript {workflow.basedir}/workflow/scripts/supplementary_CHM13_benchmarking.R \
      --rm {input.rm} \
      --cs {input.cs} \
      --hierarchy {input.hierarchy} \
      --outdir {params.outdir}
    """

BENCHMARK_SAMPLES = {
    "CHM13": ["hap1"],
    "HG002": ["hap1", "hap2"],
}

BENCHMARK_FEATURESETS = [
    "chromosome", "region", "subtelomeric",
    "acrocentric", "repeat", "gene",
]


def get_assembly_benchmark_files(wildcards):
    """All benchmark files that must exist before aggregation can run.

    Expands to (1 get_featureIDs_assembly + N featuresets * smooth_features_assembly)
    per sample/hap.
    """
    db = wildcards.database
    files = []
    for sample, haps in BENCHMARK_SAMPLES.items():
        for hap in haps:
            files.append(
                f"results/{sample}/assembly/{hap}/features/{db}/"
                f"{sample}.{hap}.{db}.benchmark.tsv"
            )
            for fs in BENCHMARK_FEATURESETS:
                files.append(
                    f"results/{sample}/assembly/{hap}/KaryoScope/{db}/"
                    f"{sample}.{hap}.{db}.{fs}.smooth_features_assembly.benchmark.tsv"
                )
    return files


rule aggregate_assembly_benchmarks:
    input:
        benchmarks=get_assembly_benchmark_files,
        script=workflow.source_path("../scripts/aggregate_assembly_benchmarks.py"),
    output:
        combined="results/benchmarks/{database}/assembly_benchmarks.combined.tsv",
        summary="results/benchmarks/{database}/assembly_benchmarks.summary.tsv",
    params:
        results_dir="results",
        samples_args=" ".join(
            f"{s}:{','.join(haps)}" for s, haps in BENCHMARK_SAMPLES.items()
        ),
        featuresets=" ".join(BENCHMARK_FEATURESETS),
    log:
        "results/benchmarks/{database}/aggregate_assembly_benchmarks.log",
    resources:
        mem_mb=5000,
        runtime=10,
    shell:
        """
        python {input.script} \
            --results-dir {params.results_dir} \
            --database {wildcards.database} \
            --samples {params.samples_args} \
            --featuresets {params.featuresets} \
            --output {output.combined} \
            --summary {output.summary} \
            2>&1 | tee {log}
        """
