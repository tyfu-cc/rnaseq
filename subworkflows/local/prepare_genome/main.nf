//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from "../../../modules/zhanglab/gunzip"
include { GUNZIP as GUNZIP_GTF              } from "../../../modules/zhanglab/gunzip"
include { GUNZIP as GUNZIP_GFF              } from "../../../modules/zhanglab/gunzip"
include { GUNZIP as GUNZIP_GENE_BED         } from "../../../modules/zhanglab/gunzip"
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from "../../../modules/zhanglab/gunzip"
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from "../../../modules/zhanglab/gunzip"

// include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
// include { UNTAR as UNTAR_RSEM_INDEX         } from '../../../modules/nf-core/untar'

include { CUSTOM_GETCHROMSIZES              } from "../../../modules/zhanglab/custom/getchromsizes"
include { GTF2BED                           } from "../../../modules/local/gtf2bed"
include { STAR_GENOMEGENERATE               } from "../../../modules/zhanglab/star/genomegenerate"
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from "../../../modules/zhanglab/rsem/preparereference"
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from "../../../modules/zhanglab/rsem/preparereference"

// include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
// include { GFFREAD                           } from '../../../modules/nf-core/gffread'
// include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'
// include { SORTMERNA as SORTMERNA_INDEX      } from '../../../modules/nf-core/sortmerna'
// include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
// include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    gene_bed                 //      file: /path/to/gene.bed
    transcript_fasta         //      file: /path/to/transcript.fasta
    star_index               // directory: /path/to/star/index/
    rsem_index               // directory: /path/to/rsem/index/
    star_rsem_index          // directory: /path/to/star/rsem/index/
    aligner                  //    string: Specifies the alignment algorithm to use
                             //            Available options: "star_rsem", :star_salmon", "hisat2"
    // gff                      //      file: /path/to/genome.gff
    // additional_fasta         //      file: /path/to/additional.fasta
    // splicesites              //      file: /path/to/splicesites.txt
    // bbsplit_fasta_list       //      file: /path/to/bbsplit_fasta_list.txt
    // sortmerna_fasta_list     //      file: /path/to/sortmerna_fasta_list.txt
    // gencode                  //   boolean: whether the genome is from GENCODE
    // featurecounts_group_type //    string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    // skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs

    main:
    ch_versions = Channel.empty()

    // Uncompress genome fasta file if required
    if (fasta.endsWith(".gz")) {
        ch_fasta = GUNZIP_FASTA ( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    }

    // Uncompress GTF annotation file or create from GFF3 if required
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith(".gz")) {
                ch_gtf = GUNZIP_GTF ( [ [:], file(gtf, checkIfExists: true) ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf, checkIfExists: true))
            }
        } else if (gff) {
            if (gff.endsWith(".gz")) {
                ch_gff      = GUNZIP_GFF ( [ [:], file(gff, checkIfExists: true) ] ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff, checkIfExists: true)).map { [ [:], it ] }
            }
            ch_gtf      = GFFREAD ( ch_gff, [] ).gtf.map { it[1] }
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }

        // Determine whether to filter the GTF or not
        // def filter_gtf =
        //     ((
        //         // Condition 1: Alignment is required and aligner is set
        //         !skip_alignment && aligner
        //     ) ||
        //     (
        //         // Condition 2: Pseudoalignment is required and pseudoaligner is set
        //         !skip_pseudo_alignment && pseudo_aligner
        //     ) ||
        //     (
        //         // Condition 3: Transcript FASTA file is not provided
        //         !transcript_fasta
        //     )) &&
        //     (
        //         // Condition 4: --skip_gtf_filter is not provided
        //         !skip_gtf_filter
        //     )
        // if (filter_gtf) {
        //     GTF_FILTER ( ch_fasta, ch_gtf )
        //     ch_gtf = GTF_FILTER.out.genome_gtf
        //     ch_versions = ch_versions.mix(GTF_FILTER.out.versions)
        // }
    }

    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    // def biotype = gencode ? "gene_type" : featurecounts_group_type
    // if (additional_fasta) {
    //     if (additional_fasta.endsWith('.gz')) {
    //         ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], file(additional_fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
    //         ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
    //     } else {
    //         ch_add_fasta = Channel.value(file(additional_fasta, checkIfExists: true))
    //     }

    //     CUSTOM_CATADDITIONALFASTA (
    //         ch_fasta.combine(ch_gtf).map { fasta, gtf -> [ [:], fasta, gtf ] },
    //         ch_add_fasta.map { [ [:], it ] },
    //         biotype
    //     )
    //     ch_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta.map { it[1] }.first()
    //     ch_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf.map { it[1] }.first()
    //     ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
    // }

    // Uncompress gene BED annotation file or create from GTF if required
    if (gene_bed) {
        if (gene_bed.endsWith(".gz")) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], file(gene_bed, checkIfExists: true) ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(gene_bed, checkIfExists: true))
        }
    } else {
        ch_gene_bed = GTF2BED( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    // Uncompress transcript fasta file / create if required
    if (transcript_fasta) {
        if (transcript_fasta.endsWith(".gz")) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], file(transcript_fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta, checkIfExists: true))
        }
        // if (gencode) {
        //     PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
        //     ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
        //     ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        // }
    } else {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA( ch_fasta, ch_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    // Create chromosome sizes file
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    // Get list of indices that need to be created
    def prepare_tool_indices = []
    if (aligner) { prepare_tool_indices << aligner }

    // Load/uncompress STAR and RSEM index or generate from scratch if required
    if ("star_rsem" in prepare_tool_indices) {
        ch_star_index = Channel.empty()
        ch_rsem_index = Channel.empty()
        ch_star_rsem_index = Channel.empty()
        
        if (star_rsem_index) {
            ch_star_rsem_index = Channel.value(file(star_rsem_index))
            ch_star_index      = ch_star_rsem_index
            ch_rsem_index      = ch_star_rsem_index
        } else if (star_index && rsem_index) {
            if (star_index.endsWith(".tar.gz")) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index))
            }
            if (rsem_index.endsWith(".tar.gz")) {
                ch_rsem_index = UNTAR_RSEM_INDEX ( [ [:], rsem_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_RSEM_INDEX.out.versions)
            } else {
                ch_rsem_index = Channel.value(file(rsem_index))
            }
            ch_star_rsem_index = ch_star_index.mix(ch_rsem_index)
        } else {
            ch_star_rsem_index = RSEM_PREPAREREFERENCE_GENOME( ch_fasta, ch_gtf ).index
            ch_star_index      = ch_star_rsem_index
            ch_versions        = ch_versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
        }
        
    }
    

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)

    star_index       = ch_star_index             // channel: path(star/index/)
    rsem_index       = ch_rsem_index             // channel: path(rsem/index/)
    star_rsem_index  = ch_star_rsem_index        // channel: path(rsem/star/index/)

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

    // splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    // bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    // rrna_fastas      = ch_rrna_fastas            // channel: path(sortmerna_fasta_list)
    // sortmerna_index  = ch_sortmerna_index        // channel: path(sortmerna/index/)
    // hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    // salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    // kallisto_index   = ch_kallisto_index         // channel: [ meta, path(kallisto/index/) ]

}
