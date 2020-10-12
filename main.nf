#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-indelible
========================================================================================
 https://github.com/ameynert/nf-indelible
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run ameynert/nf-indelible --input '/path/to/bams/*.bam' --ped samples.ped --outdir /path/to/outdir
                                       --reference hg38.fa --indelible /path/to/indelible
                                       --config /path/to/indelible/config.yml
                                       --frequency /path/to/indelible/frequency/database

    Mandatory arguments:
      --input                The input directory containing BAM files
      --outdir               The output directory where the results will be saved
      --ped                  Pedigree file relating samples in trios
      --reference            Reference genome FASTA file
      --indelible            Indelible virtualenv
      --config               Indelible configuration YAML file
      --frequency            Indelible frequency database

    Other options:
      -name                  Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (!params.input) {
    exit 1, "Input files not specified"
}

if (!params.outdir) {
    exit 1, "Output directory not specified"
}

/*if (!params.ped) {
    exit 1, "Pedigree file not specified"
}*/

if (!params.reference) {
    exit 1, "Reference genome not specified"
}

if (!params.indelible) {
    exit 1, "Indelible environment not specified"
}

if (!params.config) {
    exit 1, "Indelible configuration YAML file not specified"
}

/*if (!params.frequency) {
    exit 1, "Indelible frequency database not specified"
}*/

/*
 * Create a channel for input BAM files
 */
Channel
  .fromFilePairs( params.input, size: 1 )
  .ifEmpty { exit 1, "Cannot find any files matching ${params.input}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }
  .into { fetch_ch }


/* TODO - filter BAM files to proband only for fetch/count/score/blast/annotate steps, put parental BAM files matched with proband sample name through to de novo step  */

/*
 * STEP 1 - Fetch reads from BAM file
 */
process fetch {

  input:
  set val(name), file(alignment) from fetch_ch

  output:
  set val(name), file('*.reads'), file('*.bam') into aggregate_ch

  script:
  """
  source ${params.indelible}/venv/bin/activate
  indelible.py fetch \
    --i ${alignment} \
    --o ${name}.reads \
    --config ${params.config}
  deactivate
  """
}

/*
 * STEP 2 - Merge information across reads
 */
process aggregate {

  input:
  set val(name), file(reads), file(alignment) from aggregate_ch

  output:
  set val(name), file('*.counts') into score_ch

  script:
  """
  source ${params.indelible}/venv/bin/activate
  indelible.py aggregate \
    --i ${reads} \
    --b ${alignment} \
    --o ${name}.counts \
    --r ${params.reference} \
    --config ${params.config}
  deactivate
  """
}

/*
 * STEP 3 - Scores positions based on the read information and sequence context
 */
process score {

  input:
  set val(name), file(counts) from score_ch

  output:
  set val(name), file('*.scored') into blast_ch

  script:
  """
  source ${params.indelible}/venv/bin/activate
  indelible.py score \
    --i ${counts} \
    --o ${name}.scored \
    --config ${params.config}
  deactivate
  """
}

/*
 * STEP 4 - BLAST search of longer clipped segments (>20bp) to find matches elsewhere
 *          in the human genome and/or repeat database
 */
process blast {

  publishDir params.outdir, mode: 'copy'

  input:
  set val(name), file(scored) from blast_ch

  output:
  file(*) into annotate_ch

  script:
  """
  source ${params.indelible}/venv/bin/activate
  indelible.py blast \
    --i ${scored} \
    --config ${params.config}
  deactivate
  """
}
