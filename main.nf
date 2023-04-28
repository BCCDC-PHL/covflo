#!/usr/bin/env nextflow

//enable domain-specific-language 2
nextflow.enable.dsl=2

/**
---------------------------------------------------------------------------------
function definition
---------------------------------------------------------------------------------
*/

def helpMe() {
  log.info """

Overview:
Nextflow pipeline for generation of phylogenetic tree of SARS-CoV-2.
Replaces snakefile.

Usage:
nextflow run main.nf -profile conda [OPTIONS]

Mandatory arguments:
 --dir                          User's directory that contains input 'config' & 'data' folders.

Optional arguments:
 --seqs                         Multi-fasta file containing consensus sequences of interest [./data/sequences.fasta]
 --ref                          Reference genome used to align reads to during guided assembly [./config/Ref.gb]
 --meta                         File containing metadata for sequences under analysis [./data/metadata.csv]
 --conda_cache                  File system path where Conda env is to be stored [./covflo-main/work/]
 --drop_strains                 Excluded strains/ samples [./config/dropped_strains.txt]
 --keep_strains                 Strains that are included [./config/included_strains.txt]
 --colors                       Colors used in final auspice visualization [./config/colors.csv]
 --lat_long                     Sample latitudes and longitudes [./config/lat_longs.csv]
 --auspice                      Specifications for visualization in auspice (ex. title) [./config/auspice_config.json]
 --n_cutoff                     Maximum allowable percentage of N in a seq [0.15]
 --clip_option                  Clipkit algorithm to use [kpic-smart-gap]
 --bl_min                       Minimum branch length for RAxML [0.0000000001]
 --precision                    Precision of rounding for branch length [6]
 --length                       Length of branches to be collapsed using Gotree [0]
 --coalescent                   Timetree coalescent timescale model [const]
 --date_inference               Timetree node date estimates [marginal]
 --divergence_units             Units used to measure divergence in phylogeny refining step [mutations]
 --clock_rate                   Timetree clock rate [0.0008]
 --clock_std_dev                Timetree standard deviation of clock rate [0.0004]
 --inference                    Type of inference used in reconstructing ancestral seqs and mutations (Augur) [joint]
 --trans_probs_80               Input pairwise transmission probabilities (cov2clusters) min threshold 0.8 [./config/SARS-CoV-2_0.8_TransProbs.txt]
 --trans_probs_90               Input pairwise transmission probabilities (cov2clusters) min threshold 0.9 [./config/SARS-CoV-2_0.9_TransProbs.txt]
 --gen_clusts_80                Genomic clusters from previous build based on min trans prob 0.8 [./config/SARS-CoV-2_0.8_GenomicClusters.txt]
 --gen_clusts_90                Genomic clusters from previous build based on min trans prob 0.9 [./config/SARS-CoV-2_0.9_GenomicClusters.txt]
 --version                      Current covflo version number
 --help                         This usage statement
        """
}

def version() {
  log.info """
  covflo version: ${workflow.manifest.version}
  """
}

//displays help upon request
if (params.help) {
  helpMe()
  exit 0 //stop running
}

//version upon request
if (params.version) {
  version()
  exit 0
}

def header() {

return """
 ██████  ██████  ██    ██ ███████ ██       ██████        \\/ |    |/
██      ██    ██ ██    ██ ██      ██      ██    ██    \\/ / \\||/  /_/___/_
██      ██    ██ ██    ██ █████   ██      ██    ██     \\/   |/ \\/
██      ██    ██  ██  ██  ██      ██      ██    ██_\\__\\_\\   |  /_____/_
 ██████  ██████    ████   ██      ███████  ██████        \\  | /          /
                                                __ _-----`  |{,-----------~
                                                          \\ }{
                                                           }{{
                                                           {{}
                                                     , -=-~{ .-^- _
=========================================================================
data directory: ${params.work_dir}
output data: ${params.work_dir}results/
output tree: ${params.work_dir}auspice/
"""
}

/**
---------------------------------------------------------------------------------
program introduction
---------------------------------------------------------------------------------
*/

// this prints program header with mandatory input and output locations
log.info header()

/**
---------------------------------------------------------------------------------
process definition
---------------------------------------------------------------------------------
*/
process filtration {

tag "Filtering to exclude strains"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple file(sequences), file(meta)

output:
path("filtered.fasta")

"""
augur filter \
--sequences ${sequences} \
--metadata ${meta} \
--exclude ${params.drop_strains} \
--include ${params.keep_strains} \
--output filtered.fasta

"""
}

process percent {

tag "Remove seqs with N > cutoff using Goalign"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
file(filt_seqs)

output:
file("removedpercent.fasta")

"""
goalign \
--auto-detect clean seqs \
--cutoff ${params.n_cutoff} \
-t ${task.cpus} \
-o "removedpercent.fasta" \
-i ${filt_seqs}
"""
}

process replace {

tag "replace N sites with -"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
file(percent_seqs)

output:
path("replaced.fasta")

"""
goalign \
--auto-detect replace \
-e -s '[^ACTGactg]' \
-n '-' \
-t ${task.cpus} \
-o "replaced.fasta" \
-i ${percent_seqs}
"""
}

process clip {

tag "Clip sites that are uninformative"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
path(replace_seqs)

output:
path("informative.fasta")

"""
clipkit ${replace_seqs} \
-m ${params.clip_option} \
-l \
-o "informative.fasta"
"""
}

process dedup {

tag "Deduplicate identical sequences using Goalign"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
file(clip_seqs)

output:
tuple path("deduped.fasta"), path("names.dedup")

"""
goalign \
--auto-detect dedup \
-l "names.dedup" \
-t ${task.cpus} \
-o "deduped.fasta" \
-i ${clip_seqs}
"""
}

process compress {

tag "Compress identical sites from alignment using Goalign"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple path(deduped_seqs), path(dedup_names)

output:
tuple path("compressed.fasta"), path("weights")

"""
goalign \
--auto-detect compress \
-t ${task.cpus} \
-o "compressed.fasta" \
--weight-out "weights" \
-i ${deduped_seqs}
"""
}

process fasttree {

tag "Build fasttree without support & fastest"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple file(compress_seqs), file(compress_weights)

output:
path("fasttree.nwk")

"""
OMP_NUM_THREADS=${task.cpus} \
fasttree \
-nosupport \
-fastest \
-out "fasttree.nwk" \
-nt ${compress_seqs}
"""
}

process resolve {

tag "Resolve multifurcations in fasttree"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
file(fasttree_tree)

output:
file("resolvedtree.nwk")

"""
gotree resolve \
-t ${task.cpus} \
-i ${fasttree_tree} \
-o "resolvedtree.nwk"
"""
}

process branches {

tag "Rescale branch lengths; model GTR+I+R; min branch length 0.0000000001 using RAxML"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple file(resolve_tree), file(compress_seqs), file(compress_weights) //from compress process output

output:
tuple file("blscaled.raxml.bestTree"), file("*.log"), file("*.bestTreeCollapsed"),
file("*.bestModel"), file("*.startTree"), file("*.rba")

"""
raxml-ng \
--evaluate \
--blmin ${params.bl_min} \
--threads ${task.cpus} \
--prefix blscaled \
--force perf_threads \
--model GTR+I+R \
--tree ${resolve_tree} \
--msa ${compress_seqs} \
--site-weights ${compress_weights}
"""
}

process round {

tag "Modify branch lengths; precisions = 6 using Gotree"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple file(branches_tree), file(branches_log), file(branches_collapsed), file(branches_model),
file(branches_start), file(branches_rba)

output:
file("brlen_round.nwk")

"""
gotree brlen round \
-p ${params.precision} \
-t ${task.cpus} \
-i ${branches_tree} \
-o "brlen_round.nwk"
"""
}

process collapse {

tag "Collapse branches of length 0 using Gotree"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
file(round_tree)

output:
file("collapse_length.nwk")

"""
gotree collapse length \
-l ${params.length} \
-t ${task.cpus} \
-i ${round_tree} \
-o "collapse_length.nwk"
"""
}

process repopulate {

tag "Repopulate tree with identical sequences using Gotree"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple file(collapse_tree), path(deduped_seqs), path(dedup_names) //from dedup process output; only need names

output:
file("repopulate.nwk")

"""
gotree repopulate \
-t ${task.cpus} \
-g {dedup_names} \
-i ${collapse_tree} \
-o "repopulate.nwk"
"""
}

process order {

tag "Reorder nodes to ease comparison using Newick Utils"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
file(repopulate_tree)

output:
file("order.nwk")

"""
nw_order ${repopulate_tree} > order.nwk
"""
}

process refine {

tag "Refining phylogeny & estimating divergence rate using Augur"
publishDir "${params.work_dir}/results/", mode: 'copy'

/** 
- estimate timetree
- use ${params.coalescent} coalescent timescale
- estimate ${params.date_inference} node dates
- use {params.clock_rate} clock rate with ${params.clock_std_dev} deviation
- Does NOT filter tips more than IQDs from clock expectation (pruning)
- Gives new tree and node data to be used for translate and ancestral rules
- Roots to Wuhan-Hu-1 MN908947.3
*/

input:
tuple file(order_tree), file(percent_seqs), file(meta)

output:
tuple file("tree.nwk"), file("branch_lengths.json")

"""
augur refine \
--tree ${order_tree} \
--alignment ${percent_seqs} \
--metadata ${meta} \
--timetree \
--output-tree "tree.nwk" \
--coalescent ${params.coalescent} \
--clock-rate ${params.clock_rate} \
--clock-std-dev ${params.clock_std_dev} \
--date-confidence \
--date-inference ${params.date_inference} \
--output-node-data "branch_lengths.json" \
--divergence-units ${params.divergence_units} \
--root oldest
"""
}

process ancestral {

tag "Reconstructing ancestral sequences and mutations"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple file(refine_tree), file(refine_bls), file(filtration_seqs)

output:
path("nt_muts.json")

"""
augur ancestral \
--tree ${refine_tree} \
--alignment ${filtration_seqs} \
--output-node-data "nt_muts.json" \
--keep-ambiguous \
--keep-overhangs \
--inference ${params.inference}
"""
}

process translate {

tag "Translating to amino acids and producing individual alignments"
publishDir "${params.work_dir}/results/", mode: 'copy'

input:
tuple file(refine_tree), file(refine_bls), file(ancestral_nodes), file(reference)

output:
file("aa_muts.json")

"""
augur translate \
--tree ${refine_tree} \
--ancestral-sequences ${ancestral_nodes} \
--reference-sequence ${reference} \
--output-node-data "aa_muts.json"
"""
}

process export {
  tag "Exporting data files for auspice"
  publishDir "${params.work_dir}/auspice", mode: 'copy'

  input:
  tuple file(meta), file(refine_tree), file(refine_bls), file(ancestral_muts),\
  file(translate_muts)

  output:
  file("ncov_na.json")

      """
      augur export v2 \
          --tree ${refine_tree} \
          --metadata ${meta} \
          --node-data ${refine_bls} ${ancestral_muts} ${translate_muts} \
          --colors ${params.colors} \
          --lat-longs ${params.lat_long} \
          --minify-json \
          --auspice-config ${params.auspice} \
          --output ncov_na.json
      """
}

process clusters {

tag "Subsample SARS-CoV-2 genomic clusters at min probs of 0.8 & 0.9"
publishDir "${params.work_dir}/results", mode: 'copy'

input:
tuple file(refine_tree), file(refine_bls), file(order_tree)

output:
tuple file("SARS-CoV-2_0.8_GenomicClusters.txt"), file("SARS-CoV-2_0.8_ClustersSummary.csv"),\
file("SARS-CoV-2_0.8_TransProbs.txt"), file("SARS-CoV-2_0.9_GenomicClusters.txt"), \
file("SARS-CoV-2_0.9_ClustersSummary.csv"), file("SARS-CoV-2_0.9_TransProbs.txt")

script:
"""
Rscript ${projectDir}/bin/A3_0.8_subsam_cov2clusters_081021_141021.R "${refine_tree}" \
"${refine_bls}" \
"${order_tree}" \
"${params.trans_probs_80}" \
"${params.gen_clusts_80}"

Rscript ${projectDir}/bin/A3_0.9_subsam_cov2clusters_081021_141021.R "${refine_tree}" \
"${refine_bls}" \
"${order_tree}" \
"${params.trans_probs_90}" \
"${params.gen_clusts_90}"
"""
}

process condense {

tag "Collapse branches < 0.0000021, convert scale, & find min number of clusters with max edge length = 6"
publishDir "${params.work_dir}/results", mode: 'copy'

input:
file(order_tree)

output:
tuple file("tree_collapse_snp.nwk"), file("tc_cluster.tsv")


"""
collapse-minimum-branch-length-scale-to-snps.py ${order_tree} > tree_collapse_snp.nwk 

TreeCluster.py -i tree_collapse_snp.nwk -o tc_cluster.tsv -t 6 -m max_clade

"""
}

/**
---------------------------------------------------------------------------------
workflow
---------------------------------------------------------------------------------
*/

workflow {
  seq_ch = Channel.fromPath(params.seqs, checkIfExists:true)
  ref_ch = Channel.fromPath(params.ref, checkIfExists:true)
  meta_ch = Channel.fromPath(params.meta, checkIfExists:true)

  filtration(seq_ch.combine(meta_ch)) | percent | replace | clip | dedup | compress | fasttree | resolve
  branches(resolve.out.combine(compress.out)) | round | collapse 
  repopulate(collapse.out.combine(dedup.out)) | order
  refine(order.out.combine(percent.out.combine(meta_ch)))
  ancestral(refine.out.combine(filtration.out))
  translate(refine.out.combine(ancestral.out.combine(ref_ch)))
  export(meta_ch.combine(refine.out.combine(ancestral.out.combine(translate.out))))
  clusters(refine.out.combine(order.out))
  condense(order.out)
}

/**
---------------------------------------------------------------------------------
optional notification of completion
---------------------------------------------------------------------------------
*/

workflow.onComplete {

    def msg = """\
        
        SARS-CoV-2 Tree Generated on Almeida
        Script A2 ready to be run.
        
        Pipeline Execution Summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        Launch dir  : ${workflow.workDir}
        Data dir    : ${workflow.params.projectDir}
        Covflo vers : ${workflow.manifest.version} 
        Exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: 'jessica.caleta@bccdc.ca', subject: 'Test email', body: msg)
}

/**

process clean {
  tag "Removing Nextflow work directory?"
    shell:
        "rm -rfv work"
}
*/
