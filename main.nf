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

    nextflow run ameynert/nf-indelible --probands '/path/to/bams/*.{bam,bai}' --parents /path/to/parent/bams/dir
                                       --ped samples.ped --outdir /path/to/outdir
                                       --reference hg38.fa --indelible /path/to/indelible
                                       --config /path/to/indelible/config.yml
                                       --frequency /path/to/indelible/frequency/database

    Mandatory arguments:
      --probands             The input directory containing proband BAM files
      --parents              The input directory containing parent BAM files
      --ped                  Pedigree file relating samples in trios
      --outdir               The output directory where the results will be saved
      --reference            Reference genome FASTA file
      --indelible            Indelible virtualenv
      --config               Indelible configuration YAML file
      --frequency            Indelible frequency database

    Other options:
      -name                  Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
      --suffix               Optional suffix (e.g. "-ready") for sample names

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

if (!params.probands) {
    exit 1, "Input proband BAM files not specified"
}

if (!params.parents) {
    exit 1, "Input parent BAM file directory not specified"
}

if (!params.outdir) {
    exit 1, "Output directory not specified"
}

if (!params.ped) {
    exit 1, "Pedigree file not specified"
}

if (!params.reference) {
    exit 1, "Reference genome not specified"
}

if (!params.indelible) {
    exit 1, "Indelible environment not specified"
}

if (!params.config) {
    exit 1, "Indelible configuration YAML file not specified"
}

if (!params.frequency) {
    exit 1, "Indelible frequency database not specified"
}

/*
 * Create a channel for proband BAM files
 */
Channel
  .fromFilePairs( params.probands, size: 2 ) { file->file.name.replaceAll(/.bam|.bai$/,'') }
  .ifEmpty { exit 1, "Cannot find any files matching ${params.probands}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }
  .set { probands_ch }

/*
 * STEP 1 - Fetch reads from proband BAM file
 */
process fetch {

  input:
  set val(name), file(alignment) from probands_ch

  output:
  set val(name), file('*.reads'), file(alignment) into reads_ch

  script:
  """
  ${params.indelible}/indelible.py fetch \
    --i ${name}.bam \
    --o ${name}.reads \
    --config ${params.config}
  """
}

/*
 * STEP 2 - Merge information across reads
 */
process aggregate {

  input:
  set val(name), file(reads), file(alignment) from reads_ch

  output:
  set val(name), file('*.counts') into counts_ch

  script:
  """
  ${params.indelible}/indelible.py aggregate \
    --i ${reads} \
    --b ${name}.bam \
    --o ${name}.counts \
    --r ${params.reference} \
    --config ${params.config}
  """
}

/*
 * STEP 3 - Scores positions based on the read information and sequence context
 */
process score {

  input:
  set val(name), file(counts) from counts_ch

  output:
  set val(name), file('*.scored') into score_ch

  script:
  """
  ${params.indelible}/indelible.py score \
    --i ${counts} \
    --o ${name}.scored \
    --config ${params.config}
  """
}

/*
 * STEP 4 - BLAST search of longer clipped segments (>20bp) to find matches elsewhere
 *          in the human genome and/or repeat database
 */
process blast {

  input:
  set val(name), file(scored) from score_ch

  output:
  set val(name), file(scored), file('*.fasta*') into blast_ch

  script:
  """
  ${params.indelible}/indelible.py blast \
    --i ${scored} \
    --config ${params.config}
  """
}

/*
 * STEP 5 - Annotation with gene/exon and merging of BLAST results.
 */
process annotate {

  publishDir params.outdir, mode: 'copy',
    saveAs: { filename -> filename - ~/params.suffix/ }

  input:
  set val(name), file(scored), file(blast) from blast_ch

  output:
  set val(name), file('*.annotated') into annotated_ch

  script:
  """
  ${params.indelible}/indelible.py annotate \
    --i ${scored} \
    --o ${name}.annotated \
    --d ${params.frequency } \
    --config ${params.config}
  """
}

/*
 * STEP 6 - Read in the PED file and map proband name to 
 *          maternal and paternal BAM files
 */
