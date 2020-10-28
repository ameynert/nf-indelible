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

    nextflow run ameynert/nf-indelible --input '/path/to/bams/*.{bam,bai}' --parents /path/to/parent/bams/dir
                                       --ped samples.ped --outdir /path/to/outdir
                                       --reference hg38.fa --indelible /path/to/indelible
                                       --config /path/to/indelible/config.yml
                                       --frequency /path/to/indelible/frequency/database
                                       --denovo

    Mandatory arguments:
      --input                The input directory containing proband BAM files
      --parents              The input directory containing parent BAM files (required if --denovo set)
      --ped                  Pedigree file relating samples in trios (required if --denovo set)
      --outdir               The output directory where the results will be saved
      --reference            Reference genome FASTA file
      --indelible            Indelible virtualenv
      --config               Indelible configuration YAML file
      --frequency            Indelible frequency database

    Other options:
      -name                  Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
      --suffix               Optional suffix (e.g. "-ready") for sample names
      --denovo               Run de novo detection step

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
    exit 1, "Input proband BAM files not specified"
}

if (params.denovo && !params.parents) {
    exit 1, "Input parent BAM file directory not specified"
}

if (!params.outdir) {
    exit 1, "Output directory not specified"
}

if (params.denovo && !params.ped) {
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
 * Create value channels for input files
 */
config_ch = Channel.value(file(params.config))
reference_ch = Channel.value(file(params.reference))
frequency_ch = Channel.value(file(params.frequency))
indelible_ch = Channel.value(file(params.indelible))

if (params.denovo) {
  ped_ch = Channel.value(file(params.ped))
  parents_ch = Channel.value(file(params.parents))
}

/*
 * Create a channel for proband BAM files
 */
Channel
  .fromFilePairs( params.input, size: 2 ) { file->file.name.replaceAll(/.bam|.bai$/,'') }
  .ifEmpty { exit 1, "Cannot find any files matching ${params.input}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }
  .set { probands_ch }

/*
 * STEP 1 - Fetch reads from proband BAM file
 */
process fetch {

  input:
  set val(name), file(alignment) from probands_ch
  file indelible from indelible_ch
  file config from config_ch

  output:
  set val(name), file('*.reads'), file(alignment) into reads_ch

  script:
  name = name.replaceAll(/${params.suffix}/, '')
  """
  ${indelible}/indelible.py fetch \
    --i ${name}${params.suffix}.bam \
    --o ${name}.reads \
    --config ${config}
  """
}

/*
 * STEP 2 - Merge information across reads
 */
process aggregate {

  input:
  set val(name), file(reads), file(alignment) from reads_ch
  file indelible from indelible_ch
  file config from config_ch
  file reference from reference_ch

  output:
  set val(name), file('*.counts') into counts_ch

  script:
  """
  ${indelible}/indelible.py aggregate \
    --i ${reads} \
    --b ${name}${params.suffix}.bam \
    --o ${name}.counts \
    --r ${reference} \
    --config ${config}
  """
}

/*
 * STEP 3 - Scores positions based on the read information and sequence context
 */
process score {

  input:
  set val(name), file(counts) from counts_ch
  file indelible from indelible_ch
  file config from config_ch

  output:
  set val(name), file('*.scored') into score_ch

  script:
  """
  ${indelible}/indelible.py score \
    --i ${counts} \
    --o ${name}.scored \
    --config ${config}
  """
}

/*
 * STEP 4 - BLAST search of longer clipped segments (>20bp) to find matches elsewhere
 *          in the human genome and/or repeat database
 */
process blast {

  input:
  set val(name), file(scored) from score_ch
  file indelible from indelible_ch
  file config from config_ch

  output:
  set val(name), file(scored), file('*.fasta*') into blast_ch

  script:
  """
  ${indelible}/indelible.py blast \
    --i ${scored} \
    --config ${config}
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
  file indelible from indelible_ch
  file config from config_ch
  file frequency from frequency_ch

  output:
  set val(name), file('*.annotated') into annotated_ch,annotated_ch2

  script:
  """
  ${indelible}/indelible.py annotate \
    --i ${scored} \
    --o ${name}.annotated \
    --d ${frequency} \
    --config ${config}
  """
}

/* 
 * STEP	7 - Apply filters to annotated output:
 *     	    https://github.com/eugenegardner/indelible#recommended-filtering
 */
process	annotated_filter {
 
  publishDir  params.outdir, mode: 'copy'

  input:
  set val(name), file(annotated) from annotated_ch2

  output:
  set val(name), file('*.filtered') into annotated_filtered_ch  

  script:
  """ 
  filter_annotated.R ${annotated} ${name}.annotated.filtered
  """ 
}

if (params.denovo) {

/*
 * STEP 6 - Read in the PED file and map proband name to 
 *          maternal and paternal BAM files, then run
 *          denovo mutation event detection
 */
process denovo {

  publishDir  params.outdir, mode: 'copy'

  input:
  set val(name), file(annotated) from annotated_ch
  file indelible from indelible_ch
  file config from config_ch
  file ped from ped_ch
  file parents from parents_ch
 
  output:
  set val(name), file('*.denovo') into denovo_ch

  script:
  """
  paternal_id=`grep ${name} ${ped} | cut -f 3`
  maternal_id=`grep ${name} ${ped} | cut -f 4`

  ${indelible}/indelible.py denovo \
    --c ${annotated} \
    --m ${parents}/\${maternal_id}${params.suffix}.bam \
    --p ${parents}/\${paternal_id}${params.suffix}.bam \
    --o ${name}.denovo \
    --config ${config}
  """
}

/*
 * STEP 8 - Apply filters to de novo output:
 *          https://github.com/eugenegardner/indelible#recommended-filtering
 */
process denovo_filter {

  publishDir  params.outdir, mode: 'copy'

  input:
  set val(name), file(denovo) from denovo_ch

  output:
  set val(name), file('*.filtered') into denovo_filtered_ch

  script:
  """
  filter_denovo.R ${denovo} ${name}.denovo.filtered
  """
}

}
