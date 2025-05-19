/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ                 } from "../../modules/zhanglab/cat/fastq"
include { CUTADAPT                  } from "../../modules/zhanglab/cutadapt"
include { CHECKSTRANDEDNESS         } from "../../modules/zhanglab/checkstrandedness"
include { PRESEQ_LCEXTRAP           } from "../../modules/zhanglab/preseq/lcextrap"
include { PICARD_MARKDUPLICATES     } from "../../modules/zhanglab/picard/markduplicates"
include { RSEM_CALCULATEEXPRESSION  } from "../../modules/zhanglab/rsem/calculateexpression"
include { SUBREAD_FEATURECOUNTS     } from "../../modules/zhanglab/subread/featurecounts"

include { FASTQC                    } from "../../modules/nf-core/fastqc/main"
include { MULTIQC                   } from "../../modules/nf-core/multiqc/main"

include { MULTIQC_CUSTOM_BIOTYPE    } from "../../modules/local/multiqc_custom_biotype"

include { BAM_MARKDUPLICATES_PICARD } from "../../subworkflows/zhanglab/bam_markduplicates_picard"
include { BAM_RSEQC                 } from "../../subworkflows/zhanglab/bam_rseqc"

include { QUANTIFY_RSEM             } from "../../subworkflows/local/quantify_rsem"
include { methodsDescriptionText    } from "../../subworkflows/local/utils_rnaseq_pipeline"
include { biotypeInGtf              } from '../../subworkflows/local/utils_rnaseq_pipeline'

include { paramsSummaryMultiqc      } from "../../subworkflows/nf-core/utils_nfcore_pipeline"
include { softwareVersionsToYAML    } from "../../subworkflows/nf-core/utils_nfcore_pipeline"




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_biotypes_header_multiqc = file("${projectDir}/assets/multiqc/biotypes_header.txt", checkIfExists: true)

workflow RNASEQ {

    take:
    ch_samplesheet      // channel: samplesheet read in from --input
    ch_fasta            // channel: path(genome.fasta)
    ch_gtf              // channel: path(genome.gtf)
    ch_fai              // channel: path(genome.fai)
    ch_chrom_sizes      // channel: path(genome.sizes)
    ch_gene_bed         // channel: path(gene.bed)
    ch_transcript_fasta // channel: path(transcript.fasta)
    ch_star_index       // channel: path(star/index/)
    ch_rsem_index       // channel: path(rsem/index/)
    ch_star_rsem_index  // channel: path(star/rsem/index/)
    ch_versions         // channel: [ path(versions.yml) ]
    

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Prepare for merging fastq files from the same sample
    ch_samplesheet
        .branch {
            meta, fastqs ->
                single: fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }

    // MODULE: Concatenate fastq files from same sample if required
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_merged_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    // MODULE: Trim the adapter
    CUTADAPT(
        ch_merged_fastq
    )
    ch_trimmed_fastq = CUTADAPT.out.reads

    // MODULE: Set the strandedness for each sample
    CHECKSTRANDEDNESS(
        ch_trimmed_fastq,
        ch_gtf,
        ch_transcript_fasta
    )
    
    CHECKSTRANDEDNESS.out.txt
        .join(ch_trimmed_fastq)
        .map {
            meta, txt, reads ->
                def strandedness = txt.readLines().first()
                return [ meta + [ strandedness: strandedness ], reads ]
        }
        .set { ch_strand_inferred_fastq }

    // MODULE: Run FastQC on trimmed reads
    FASTQC (
        ch_strand_inferred_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // SUBWORKFLOW: Run RSEM (Do the alignment using the --star option)
    QUANTIFY_RSEM(
        ch_strand_inferred_fastq,
        ch_fasta.map { [ [:], it ] },
        ch_star_rsem_index
    )
    ch_genome_bam       = QUANTIFY_RSEM.out.bam
    ch_genome_bam_index = QUANTIFY_RSEM.out.bai
    ch_star_log         = QUANTIFY_RSEM.out.logs
    ch_multiqc_files    = ch_multiqc_files.mix(ch_star_log.collect{it[1]})

    // MODULE: Run Preseq
    if (!params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect{it[1]})
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    // SUBWORKFLOW: Run Picard MarkDuplicates
    if (!params.skip_markduplicates) {
        BAM_MARKDUPLICATES_PICARD(
            ch_genome_bam,
            ch_fasta.map { [ [:], it ] },
            ch_fai.map { [ [:], it ] }
        )
        ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_multiqc_files    = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]})
        // ch_multiqc_files    = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]})
        // ch_multiqc_files    = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]})
        // ch_multiqc_files    = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]})

        //if (params.bam_csi_index) {
        //    ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
        //}
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }

    // MODULES: Run featureCounts
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    if (!params.skip_biotype_qc && biotype) {

        ch_gtf
            .map { biotypeInGtf(it, biotype) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .combine(ch_gtf)
            .combine(biotype_in_gtf)
            .filter { it[-1] }
            .map { it[0..<it.size()-1] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header_multiqc
        )
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv.collect{it[1]})
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }

    // Get RSeqC modules to run
    def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(",").collect{ it.trim().toLowerCase() } : []
    if (params.bam_csi_index) {
        for (rseqc_module in ["read_distribution", "inner_distance", "tin"]) {
            if (rseqc_modules.contains(rseqc_module)) {
                rseqc_modules.remove(rseqc_module)
            }
        }
    }
    if (!params.skip_rseqc && rseqc_modules.size() > 0) {

        BAM_RSEQC (
            ch_genome_bam.join(ch_genome_bam_index, by: [0]),
            ch_gene_bed,
            rseqc_modules
        )

        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.genebody_coverage_txt.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.genebody_coverage_pdf.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.genebody_coverage_rscript.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.innerdistance_freq.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionannotation_log.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readdistribution_txt.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readgc_xls.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readgc_rscript.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readgc_pdf.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.fragsize_txt.collect{it[1]})
        ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)
    }

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'rnaseq-nf_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // MODULE: MultiQC
    ch_multiqc_config        = Channel.fromPath(
        "${projectDir}/assets/multiqc/multiqc_config.yaml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )


    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
