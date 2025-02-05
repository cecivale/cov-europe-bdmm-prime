# Modified Nextrain workflow from 2020-09-11 by Cecilia Valenzuela
# CBB Master's Thesis A comprehensive study of the earllly phylodynamics of SARS-CoV-2 in Europe
# cEvo group DBSSE ETH 
# MIT License
# Copyright (c) 2020 Nextstrain, full text https://github.com/nextstrain/ncov/blob/master/LICENSE

import pandas as pd
import random
from snakemake.utils import validate
from datetime import date
from treetime.utils import numeric_date

from itertools import product, chain

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

BUILD_NAMES = list(config["builds"].keys())
ANALYSIS_NAME = list(chain.from_iterable([config["builds"][build_name]["beast_analysis"] for build_name in BUILD_NAMES]))
SEED = 1234
random.seed(SEED)

# Define patterns we expect for wildcards.
wildcard_constraints:
    build_name = r"[a-zA-Z0-9_.-]+",
    analysis_name = r"[a-zA-Z0-9_.-]+",
    date = r"[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]",
    seed = r"[0-9]*"

localrules: load

def filter_combinator(combinator, combinations):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            # Use frozenset instead of tuple
            # in order to accomodate
            # unpredictable wildcard order
            if frozenset(wc_comb) in combinations:
                yield wc_comb
    return filtered_combinator

combinations = list(chain.from_iterable([[frozenset(sett) for sett in product(frozenset({("build_name", build_name)}), frozenset(product(["analysis_name"], config["builds"][build_name]["beast_analysis"])))] for build_name in BUILD_NAMES]))
filtered_product = filter_combinator(product, combinations)
PARTICLES = config["beast"]["n_particles"]

rule all:
    input:
        ml_tree = expand("results/{build_name}/tree.nwk", build_name=BUILD_NAMES),
        summary_tree = expand(expand("results/{build_name}/{analysis_name}.{{particles}}.mcc.typed.node.tree", filtered_product, build_name=BUILD_NAMES, analysis_name=ANALYSIS_NAME), particles=PARTICLES) ,
        report_traj = expand(expand("results/{build_name}/{analysis_name}.{{particles}}.analyze_trajectories.pdf", filtered_product, build_name=BUILD_NAMES, analysis_name=ANALYSIS_NAME), particles=PARTICLES),
        #figs_traj = expand(expand("results/{build_name}/figs_{analysis_name}_{{particles}}", filtered_product, build_name=BUILD_NAMES, analysis_name=ANALYSIS_NAME), particles=PARTICLES),
        #tables_traj = expand("results/{build_name}/tables_{analysis_name}", filtered_product, build_name=BUILD_NAMES, analysis_name=ANALYSIS_NAME),
        summary_log =  expand("results/{build_name}/{analysis_name}_comb.summary.tsv", filtered_product, build_name=BUILD_NAMES, analysis_name=ANALYSIS_NAME),
        

rule clean:
    message: "Removing directories: {params}"
    params:
        "results",
        "logs",
        "benchmarks"
    shell:
        "rm -rfv {params}"

# Include small, shared functions that help build inputs and parameters.
include: "common.smk"


rule load:
    message: "Loading sequences and metadata from files"
    output:
        sequences = config["sequences"],
        metadata = config["metadata"]


rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.load.output.sequences,
        metadata = rules.load.output.metadata,
        include = config["files"]["include"],
        exclude = config["files"]["exclude"]
    output:
        sequences = "results/filtered.fasta"
    log:
        "logs/filtered.txt"
    params:
        min_length = config["filter"]["min_length"],
        exclude_where = config["filter"]["exclude_where"],
        min_date = config["filter"]["min_date"],
        date = config["filter"]["max_date"]
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --max-date {params.date} \
            --min-date {params.min_date} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where}\
            --min-length {params.min_length} \
            --output {output.sequences} 2>&1 | tee {log}
        """


# I removed rule excluded_sequences, rule align_excluded, rule diagnose_excluded from 
# nextstrain workflow, since I think the information on log file about excluded sequences 
# from rule filter is enough for our analysis.


checkpoint partition_sequences:
    input:
        sequences = rules.filter.output.sequences
    output:
        split_sequences = directory("results/split_sequences/")
    log:
        "logs/partition_sequences.txt"
    params:
        sequences_per_group = config["partition_sequences"]["sequences_per_group"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/partition-sequences.py \
            --sequences {input.sequences} \
            --sequences-per-group {params.sequences_per_group} \
            --output-dir {output.split_sequences} 2>&1 | tee {log}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - gaps relative to reference are considered real
        Cluster:  {wildcards.cluster}
        """
    input:
        sequences = "results/split_sequences/{cluster}.fasta",
        reference = config["files"]["reference"]
    output:
        alignment = "results/split_alignments/{cluster}.fasta"
    log:
        "logs/align_{cluster}.txt"
    benchmark:
        "benchmarks/align_{cluster}.txt"
    threads: 2
    conda: config["conda_environment"]
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference 2>&1 | tee {log}
        """

def _get_alignments(wildcards):
    checkpoint_output = checkpoints.partition_sequences.get(**wildcards).output[0]
    return expand("results/split_alignments/{i}.fasta",
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)

rule aggregate_alignments:
    message: "Collecting alignments"
    input:
        alignments = _get_alignments
    output:
        alignment = "results/aligned.fasta"
    log:
        "logs/aggregate_alignments.txt"
    conda: config["conda_environment"]
    shell:
        """
        cat {input.alignments} > {output.alignment} 2> {log}
        """

rule diagnostic:
    message: "Scanning aligned sequences {input.alignment} for problematic sequences"
    input:
        alignment = rules.aggregate_alignments.output.alignment,
        metadata = rules.load.output.metadata,
        reference = config["files"]["reference"]
    output:
        diagnostics = "results/sequence-diagnostics.tsv",
        flagged = "results/flagged-sequences.tsv",
        to_exclude = "results/to-exclude.txt"
    log:
        "logs/diagnostics.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/diagnostic.py \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --reference {input.reference} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --output-flagged {output.flagged} \
            --output-diagnostics {output.diagnostics} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """


rule refilter:
    message:
        """
        excluding sequences flagged in the diagnostic step in file {input.exclude}
        """
    input:
        sequences = rules.aggregate_alignments.output.alignment,
        metadata = rules.load.output.metadata,
        exclude = rules.diagnostic.output.to_exclude
    output:
        sequences = "results/aligned-filtered.fasta"
    log:
        "logs/refiltered.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} 2>&1 | tee {log}
        """


rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.refilter.output.sequences
    output:
        alignment = "results/masked.fasta"
    log:
        "logs/mask.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"],
        mask_sites = config["mask"]["mask_sites"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --mask-terminal-gaps \
            --output {output.alignment} 2>&1 | tee {log}
        """

def _get_subsampling_settings(wildcards):
    # Allow users to override default subsampling with their own settings keyed
    # by location type and name. For example, "region_europe" or
    # "country_iceland". Otherwise, default to settings for the location type.
    subsampling_scheme = _get_subsampling_scheme_by_build_name(wildcards.build_name)
    subsampling_settings = config["subsampling"][subsampling_scheme]

    if hasattr(wildcards, "subsample"):
        return subsampling_settings[wildcards.subsample]
    else:
        return subsampling_settings


def get_priorities(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return f"results/{wildcards.build_name}/proximity_{subsampling_settings['priorities']['focus']}.tsv"
    else:
        # TODO: find a way to make the list of input files depend on config
        return config["files"]["include"]


def get_priority_argument(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return "--priority " + get_priorities(wildcards)
    else:
        return ""


def _get_specific_subsampling_setting(setting, optional=False):
    def _get_setting(wildcards):
        if optional:
            value = _get_subsampling_settings(wildcards).get(setting, "NA")
        else:
            value = _get_subsampling_settings(wildcards)[setting]

        if isinstance(value, str):
            # Load build attributes including geographic details about the
            # build's region, country, division, etc. as needed for subsampling.
            build = config["builds"][wildcards.build_name]
            value = value.format(**build)
        else:
            return value

        # Check format strings that haven't been resolved.
        if re.search(r'\{.+\}', value):
            raise Exception(f"The parameters for the subsampling scheme '{wildcards.subsample}' of build '{wildcards.build_name}' reference build attributes that are not defined in the configuration file: '{value}'. Add these build attributes to the appropriate configuration file and try again.")

        return value

    return _get_setting


# rule subsample:
#     message:
#         """
#         Subsample all sequences into a {wildcards.subsample} set for build '{wildcards.build_name}' with {params.sequences_per_group} sequences.
#         """
#     input:
#         sequences = rules.mask.output.alignment,
#         metadata = rules.load.output.metadata,
#         include = config["files"]["include"],
#         priorities = get_priorities
#     output:
#         sequences = "results/{build_name}/sample-{subsample}.fasta"
#     params:
#         group_by = _get_specific_subsampling_setting("group_by"),
#         sequences_per_group = _get_specific_subsampling_setting("seq_per_group"),
#         #max_sequences = _get_specific_subsampling_setting("max_sequences", optional=True),
#         exclude_argument = _get_specific_subsampling_setting("exclude", optional=True),
#         include_argument = _get_specific_subsampling_setting("include", optional=True),
#         query_argument = _get_specific_subsampling_setting("query", optional=True),
#         priority_argument = get_priority_argument
#     log:
#         "logs/subsample_{build_name}_{subsample}.txt"
#     conda: config["conda_environment"]
#     shell:
#         """
#         augur filter \
#             --sequences {input.sequences} \
#             --metadata {input.metadata} \
#             --max-date 2020-03-08 \
#             --include {input.include} \
#             {params.exclude_argument} \
#             {params.include_argument} \
#             {params.query_argument} \
#             {params.priority_argument} \
#             --group-by {params.group_by} \
#             --sequences-per-group {params.sequences_per_group} \
#             --subsample-seed SEED \
#             --output {output.sequences} 2>&1 | tee {log}
#         """

rule rsubsample:
    message:
        """
        R weighted subsampling all sequences into a {wildcards.subsample} set for build '{wildcards.build_name}' with {params.seq_per_deme} sequences.
        """
    input:
        sequences = rules.mask.output.alignment,
        metadata = rules.load.output.metadata,
        include = config["files"]["include"]
    output:
        sequences = "results/{build_name}/sample-{subsample}.fasta",
        fig = "results/{build_name}/sample-figs/sample-{subsample}.png"
    params:
        region = _get_specific_subsampling_setting("region", optional=True),
        country = _get_specific_subsampling_setting("country", optional=True),
        division = _get_specific_subsampling_setting("division", optional=True),
        exclude_country = _get_specific_subsampling_setting("exclude_country", optional=True),
        exclude_division = _get_specific_subsampling_setting("exclude_division", optional=True),
        min_date = _get_specific_subsampling_setting("min_date", optional=False),
        max_date = _get_specific_subsampling_setting("max_date", optional=False),
        seq_per_deme = _get_specific_subsampling_setting("seq_per_deme", optional=False),
        seed = sum([ord(v) for i,v in enumerate("{build_name}")]) # Using build name to set a seed
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    conda: config["conda_environment"]
    shell:
        """
        Rscript scripts/subsample.R \
            --alignment {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --region {params.region} \
            --country {params.country} \
            --division {params.division} \
            --exclude_country {params.exclude_country} \
            --exclude_division {params.exclude_division} \
            --from {params.min_date} \
            --to {params.max_date} \
            --seq_per_deme {params.seq_per_deme} \
            --seed {params.seed} \
            --prob "cases" \
            --output {output.sequences} \
            --output_figure {output.fig} 2>&1 | tee {log}
        """
        

# Not sure if we are interested in this
rule proximity_score:
    message:
        """
        determine priority for inclusion in as phylogenetic context by
        genetic similiarity to sequences in focal set for build '{wildcards.build_name}'.
        """
    input:
        alignment = rules.mask.output.alignment,
        metadata = rules.load.output.metadata,
        reference = config["files"]["reference"],
        focal_alignment = "results/{build_name}/sample-{focus}.fasta"
    output:
        priorities = "results/{build_name}/proximity_{focus}.tsv"
    log:
        "logs/subsampling_priorities_{build_name}_{focus}.txt"
    resources: mem_mb = 4000
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/priorities.py --alignment {input.alignment} \
            --metadata {input.metadata} \
            --reference {input.reference} \
            --focal-alignment {input.focal_alignment} \
            --output {output.priorities} 2>&1 | tee {log}
        """

def _get_subsampled_files(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    return [
        f"results/{wildcards.build_name}/sample-{subsample}.fasta"
        for subsample in subsampling_settings
    ]

rule combine_samples:
    message:
        """
        Combine and deduplicate FASTAs
        """
    input:
        _get_subsampled_files
    output:
        alignment = "results/{build_name}/subsampled_alignment.fasta"
    log:
        "logs/subsample_regions_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py \
            --input {input} \
            --output {output} 2>&1 | tee {log}
        """

rule adjust_names:
    message: 
        """
        Adjust metadata and final alignment sequence names
        """
    input: 
        alignment = rules.combine_samples.output.alignment,
        metadata = rules.load.output.metadata,
    output: 
        alignment = "results/{build_name}/final_alignment.fasta",
        metadata = "results/{build_name}/final_metadata.tsv"
    params:
        demes = _get_subsampling_settings
    log:
        "logs/{build_name}_adjust_names.txt"
    shell:
        """
        Rscript scripts/adjust_names.R \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --demes {params.demes} \
            --output_alignment {output.alignment} \
            --output_metadata {output.metadata} 2>&1 | tee {log}
        """


# Build a preliminary tree
rule tree:
    message: 
        """
        Building preliminary tree with augur
        """
    input:
        alignment = rules.adjust_names.output.alignment
    output:
        tree = "results/{build_name}/tree_raw.nwk"
    params:
        args = lambda w: config["tree"].get("tree-builder-args","") if "tree" in config else ""
    log:
        "logs/tree_{build_name}.txt"
    benchmark:
        "benchmarks/tree_{build_name}.txt"
    threads: 16
    resources:
        # Multiple sequence alignments can use up to 40 times their disk size in
        # memory, especially for larger alignments.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 40 * int(input.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            --nthreads {threads} 2>&1 | tee {log}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.adjust_names.output.alignment,
        #metadata = _get_metadata_by_wildcards
        metadata = rules.adjust_names.output.metadata
    output:
        tree = "results/{build_name}/tree.nwk",
        node_data = "results/{build_name}/branch_lengths.json"
    log:
        "logs/refine_{build_name}.txt"
    benchmark:
        "benchmarks/refine_{build_name}.txt"
    threads: 1
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    params:
        root = config["refine"]["root"],
        clock_rate = config["refine"]["clock_rate"],
        clock_std_dev = config["refine"]["clock_std_dev"],
        coalescent = config["refine"]["coalescent"],
        date_inference = config["refine"]["date_inference"],
        divergence_unit = config["refine"]["divergence_unit"],
        clock_filter_iqd = config["refine"]["clock_filter_iqd"],
        keep_polytomies = "--keep-polytomies" if config["refine"].get("keep_polytomies", False) else "",
        timetree = "" if config["refine"].get("no_timetree", False) else "--timetree"
    conda: config["conda_environment"]
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            {params.timetree} \
            {params.keep_polytomies} \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --divergence-unit {params.divergence_unit} \
            --date-confidence \
            --no-covariance \
            --clock-filter-iqd {params.clock_filter_iqd} 2>&1 | tee {log}
        """


rule beast:
    message: 
        """
        Running BEAST2 analysis BDMM-Prime, MCMC chain {wildcards.seed}.
        """
    input:
        #xml = _get_beast_analysis,
        xml = "analyses/{analysis_name}.xml",
        alignment = rules.adjust_names.output.alignment
    output:
        trace = "results/{build_name}/chains/{analysis_name}.{seed}_f.log",
        trees = "results/{build_name}/chains/{analysis_name}.{seed}_f.trees",
    params:
        log = "beast_{build_name}_{analysis_name}.{seed}.txt",
        jar = config["beast"]["jar"],
        length = config["beast"]["l_mcmc"],
        logevery = round(config["beast"]["l_mcmc"]/10000),
        seed = config["beast"]["n_mcmc"],
        name = "results/{build_name}/chains/{analysis_name}.{seed}",
        action = "overwrite"
    benchmark:
        "benchmarks/beast_{build_name}_{analysis_name}.{seed}.benchmark.txt"
    resources:
        runtime = config["beast"]["t_mcmc"],
        mem_mb = 4096,
        cpus = config["beast"]["cpus_mcmc"]
    shell:
        """
        scp {input.xml} results/{wildcards.build_name}/{wildcards.analysis_name}.xml
        cd results/{wildcards.build_name}/ 
        mkdir -p chains
        java -Xmx3G -jar {params.jar} -D "chain_length={params.length}" -D "log_every={params.logevery}" -seed {wildcards.seed} -statefile {wildcards.analysis_name}.{wildcards.seed}.state -{params.action} {wildcards.analysis_name}.xml 2>&1 | tee -a {params.log} || :
        cd ../../
        mv results/{wildcards.build_name}/{params.log} logs/

        touch {params.name}.log {params.name}.trees
        scp {params.name}.log {output.trace}
        scp {params.name}.trees {output.trees}
        #find . -type f -empty -delete
        """


rule summarize_trace:
    message: 
        """
        Summarizing trace file {input.trace} with Log Analyser v1.8.2.
        """
    input: 
        trace = rules.beast.output.trace
    output:
        trace_summary = "results/{build_name}/chains/{analysis_name}.{seed}.summary.txt"
    log:
        "logs/summarize_trace_{build_name}_{analysis_name}.{seed}.txt"
    params:
        jar = config["beast"]["jar"],
        burnin = config["beast"]["burnin"]
    shell:
        """
        java -cp {params.jar} beast.util.LogAnalyser -b {params.burnin} {input.trace} > {output.trace_summary} 2> {log} 
        # java -cp {params.jar} beast.util.LogAnalyser -b {params.burnin} {input.trace} > {output.trace_summary} 2> {log} || :
        # touch {output.trace_summary}
        """

def _get_tracesummary(wildcards):
    files = expand("results/{{build_name}}/chains/{{analysis_name}}.{seed}.summary.txt",
        seed=config["beast"]["n_mcmc"])
    return files

checkpoint chains_diagnostic:
    message: 
        """
        Check that trace file have all ESS values >= {params.min_ess} to include it in the combined chain.
        """
    input:
        trace_summary = _get_tracesummary
    output:
        diagnostic = "results/{build_name}/{analysis_name}_chains_diagnostic.txt",
    params:
        min_ess = config["beast"]["min_ess"],
        burnin = config["beast"]["burnin"],
        chain_length = config["beast"]["l_mcmc"]
    log:
        "logs/chains_diagnostic_{build_name}_{analysis_name}.txt"
    shell:
        """
        Rscript scripts/chains_diagnostic.R \
        --input {input.trace_summary} \
        --ess {params.min_ess} \
        --burnin {params.burnin} \
        --length {params.chain_length} \
        --output {output.diagnostic} 2>&1 | tee -a {log}
        """

def _get_trace_tocombine(wildcards):
    diagnostic = pd.read_csv(checkpoints.chains_diagnostic.get(**wildcards).output.diagnostic, sep="\t")
    files = expand(
        expand("results/{build_name}/chains/{analysis_name}.{{seed}}_f.log", zip,
        build_name=wildcards.build_name, analysis_name=wildcards.analysis_name),
        seed=diagnostic[diagnostic["included"] == 1]["seed"])
    if not files:
        files = "message_trace_" + wildcards.build_name + "_" + wildcards.analysis_name + ".txt"
        with open(files, "w") as file:
            file.write("None of the trace files have the indicated minimum of ESS for all variables. Change min_ESS in config file if you want to combine anyway, or resume beast chains by running snakemake -R beast --config action='resume'")
    return files


def _get_trees_tocombine(wildcards):
    diagnostic = pd.read_csv(checkpoints.chains_diagnostic.get(**wildcards).output.diagnostic, sep="\t")
    files = expand(
        expand("results/{build_name}/chains/{analysis_name}.{{seed}}_f.trees", zip,
        build_name=wildcards.build_name, analysis_name=wildcards.analysis_name),
        seed=diagnostic[diagnostic["included"] == 1]["seed"])
    if not files:
        files = "message_trees_" + wildcards.build_name + "_" + wildcards.analysis_name + ".txt"
        with open(files, "w") as file:
            file.write("None of the trace files have the indicated minimum of ESS for all variables. Change min_ESS in config file if you want to combine anyway, or resume beast chains by running snakemake -R beast")
    return files


rule combine_trace:
    message: 
        """
        Combine trace files with ESS >= {params.min_ess}: {input.trace_files} with LogCombiner v1.8.2.
        """
    input:
        trace_files = _get_trace_tocombine    
    output:
        combined_trace = "results/{build_name}/{analysis_name}_comb.log"
    log:
        "logs/combine_trace_{build_name}_{analysis_name}.txt"
    params:
        jar = config["beast"]["jar"],
        burnin = config["beast"]["burnin"],
        min_ess = config["beast"]["min_ess"],
        empty_message = "message_trace_{build_name}_{analysis_name}.txt",
        input_command = lambda wildcards, input: " -log ".join(input) 
        #input_command = lambda wildcards, input: [" -log  " + f for f in input]
    shell:
        """
        if [ -e {params.empty_message} ] 
        then 
            cat {params.empty_message} | tee {log} && rm {params.empty_message}
        else
            java -cp {params.jar} beast.app.tools.LogCombinerLauncher -log {params.input_command} -o {output.combined_trace} -b {params.burnin}  2>&1 | tee -a {log}
        fi
        """

rule combine_trees:
    message: 
        """
        Combine tree files whose trace files with ESS >= {params.min_ess}: {input.tree_files} with LogCombiner v1.8.2
        """
    input:
        tree_files = _get_trees_tocombine
    output:
        combined_trees = "results/{build_name}/{analysis_name}_comb.trees"
    log:
        "logs/combine_trees_{build_name}_{analysis_name}.txt"
    params:
        jar = config["beast"]["jar"],
        burnin = config["beast"]["burnin"],
        min_ess = config["beast"]["min_ess"],
        empty_message = "message_trees_{build_name}_{analysis_name}.txt",
        input_command = lambda wildcards, input: " -log ".join(input) 
    shell:
        """
        if [ -e {params.empty_message} ] 
        then 
            cat {params.empty_message} | tee {log} && rm {params.empty_message}
        else
            java -cp {params.jar} beast.app.tools.LogCombinerLauncher -log {params.input_command} -o {output.combined_trees}  -b {params.burnin} 2>&1 | tee -a {log}
        fi
        """

rule resample_trace:
    input:
        trace = rules.combine_trace.output.combined_trace
    output:
        thinned_trace = "results/{build_name}/{analysis_name}_comb.thinned.log"
    log:
        "logs/resample_trace_{build_name}_{analysis_name}.txt"
    params:
        jar = config["beast"]["jar"],
        resample = config["beast"]["resample"]
    shell:
        """
        java -cp {params.jar} beast.app.tools.LogCombinerLauncher -log {input.trace} -o {output.thinned_trace} -resample {params.resample} 2>&1 | tee {log}
        """

rule resample_trees:
    input:
        trees = rules.combine_trees.output.combined_trees
    output:
        thinned_trees = "results/{build_name}/{analysis_name}_comb.thinned.trees"
    log:
        "logs/resample_trees_{build_name}_{analysis_name}.txt"
    params:
        jar = config["beast"]["jar"],
        resample = config["beast"]["resample"]
    shell:
        """
        java -cp {params.jar} beast.app.tools.LogCombinerLauncher -log {input.trees} -o {output.thinned_trees} -resample {params.resample}  2>&1 | tee {log}
        """

rule trajectory_mapping:
    input: 
        xml = "analyses/{analysis_name}.trajectoryMapper.xml",
        trace = rules.resample_trace.output.thinned_trace,
        trees = rules.resample_trees.output.thinned_trees
    output:
        typed_trees = "results/{build_name}/{analysis_name}.{particles}.typed.trees",
        node_typed_trees = "results/{build_name}/{analysis_name}.{particles}.typed.node.trees",
        trajectories = "results/{build_name}/{analysis_name}.{particles}.TL.traj"
    params:
        log = "trajectory_mapping_{build_name}_{analysis_name}_{particles}.txt",
        jar = config["beast"]["jar"],
        particles = config["beast"]["n_particles"]
    benchmark:
        "benchmarks/trajectory_mapping_{build_name}_{analysis_name}_{particles}.benchmark.txt"
    resources:
        runtime = config["beast"]["t_traj"],
        mem_mb = 4096,
        cpus = 2
    shell:
        """
        scp {input.xml} results/{wildcards.build_name}/{wildcards.analysis_name}.xml
        cd results/{wildcards.build_name}/
        java -Xmx3G -jar {params.jar} -D "nParticles={wildcards.particles}" -overwrite {wildcards.analysis_name}.xml 2>&1 | tee {params.log}  || :
        cd ../../
        mv results/{wildcards.build_name}/{params.log} logs/
        touch {output.typed_trees} {output.node_typed_trees} {output.trajectories}
        """

rule clean_beast:
    input:
        xml = "results/{build_name}/{analysis_name}.xml",
        trajectories = rules.trajectory_mapping.output.trajectories
    shell:
        """
        rm results/{wildcards.build_name}/{wildcards.analysis_name}.xml 
        find results/{wildcards.build_name}/ -type f -empty -delete
        """

rule summarize_comb_trace:
    message: 
        """
        Create summary table for combined trace {input.trace} with LogAnalyser v1.8.2.
        """
    input:
        trace = rules.combine_trace.output.combined_trace
    output:
        summary_trace = "results/{build_name}/{analysis_name}_comb.summary.tsv"
    log:
        "logs/summarize_comb_trace_{build_name}_{analysis_name}.txt"
    params:
        jar = config["beast"]["jar"]
    shell:
        """
        java -cp {params.jar} beast.util.LogAnalyser {input.trace} > {output.summary_trace} 2> {log}
        """


rule summarize_trees:
    message: 
        """
        Summarize trees to Maximum clade credibility tree with mean node heights with TreeAnnotator v1.8.2.
        """
    input:
        combined_trees = rules.trajectory_mapping.output.node_typed_trees
    output:
        summary_tree = "results/{build_name}/{analysis_name}.{particles}.mcc.typed.node.tree"
    log:
        "logs/summarize_trees_{build_name}_{analysis_name}_{particles}.txt"
    params:
        jar = config["beast"]["jar"],
        heights = config["beast"]["mcc_heights"]
    shell:
        """
        java -cp {params.jar} beast.app.treeannotator.TreeAnnotatorLauncher -heights {params.heights} {input.combined_trees} {output.summary_tree} 2>&1 | tee {log} 
        """


rule analyze_trajectories:
    message: 
        """
        Analyze and plot trajectories from {input.trajectories}.
        """
    input:
        trajectories = rules.trajectory_mapping.output.trajectories,
        metadata = rules.adjust_names.output.metadata
    output:
        report = "results/{build_name}/{analysis_name}.{particles}.analyze_trajectories.pdf"
    log:
        "logs/analyze_trajectories_{build_name}_{analysis_name}_{particles}.txt"
    benchmark:
        "benchmarks/analyze_trajectories_{build_name}_{analysis_name}_{particles}.benchmark.txt"
    params:
        burnin = 0, #config["beast"]["burnin"],
        demes = _get_subsampling_settings,
        n_traj =  config["beast"]["n_traj"],
        figs = "results/{build_name}/figs_{analysis_name}_{particles}",
        tables = "results/{build_name}/tables_{analysis_name}_{particles}"
    shell:
        """
        mkdir -p {params.figs} {params.tables}
        Rscript -e "rmarkdown::render('scripts/analyze_trajectories.Rmd', output_file='../{output.report}')" \
            --args = {input.trajectories} \
            --burnin = {params.burnin} \
            --metadata = {input.metadata} \
            --demes = {params.demes} \
            --n = {params.n_traj} \
            --output_figure = {params.figs}/{wildcards.analysis_name}_figtraj \
            --output_table = {params.tables}/{wildcards.analysis_name}_table 2>&1 | tee -a {log}
        """


# rule get_sampling_proportions:
#     params:
#         demes = str(config["subsampling"]["ssprop_China"])
#     output:
#         sampling_proportions = "results/sampling_proportions_ssprop_China.tsv"
#     shell:
#         """
#         Rscript scripts/sampling_proportion.R --input {params.demes} --output {output.sampling_proportions}
#         """
