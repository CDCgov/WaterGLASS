
include { BBMAP_ALIGN_TO_REF           } from '../../modules/nf-core/bbmap/align_withRef/main'
include { BBMAP_ALIGN_TO_SKESA         } from '../../modules/nf-core/bbmap/align_toskesa/main'
include { BBMAP_ALIGN_TO_UNICYCLER     } from '../../modules/nf-core/bbmap/align_tounicycler/main'

workflow BBMAP_ALIGN{

    take:
    ch_assembly_fastq         // trimmed fastq reads
    fasta                   //  reference genome
    ch_assemblies       // unicycler assemblies
    ch_skesa_assemblies          // skesa assemblies



    main:
    ch_versions = Channel.empty()
    ch_unicycler_bbmap = Channel.empty()
    ch_skesa_bbmap = Channel.empty()

    ch_unicycler_bbmap  = ch_assembly_fastq.join(ch_assemblies)
    ch_skesa_bbmap      = ch_assembly_fastq.join(ch_skesa_assemblies)

    BBMAP_ALIGN_TO_REF( ch_assembly_fastq, fasta )
 
    ch_versions = ch_versions.mix(BBMAP_ALIGN_TO_REF.out.versions)


    BBMAP_ALIGN_TO_SKESA(ch_skesa_bbmap)

    ch_versions = ch_versions.mix(BBMAP_ALIGN_TO_SKESA.out.versions)


    BBMAP_ALIGN_TO_UNICYCLER(ch_unicycler_bbmap)

    ch_versions = ch_versions.mix(BBMAP_ALIGN_TO_UNICYCLER.out.versions)

    emit:

    bam_ref         =   BBMAP_ALIGN_TO_REF.out.bam
    bam_skesa       =   BBMAP_ALIGN_TO_SKESA.out.bam
    bam_unicycler   =   BBMAP_ALIGN_TO_UNICYCLER.out.bam

    ref_log         =   BBMAP_ALIGN_TO_REF.out.log.collect { it[1] }
    skesa_log       =   BBMAP_ALIGN_TO_SKESA.out.log.collect { it[1] }
    unicycler_log   =   BBMAP_ALIGN_TO_UNICYCLER.out.log.collect { it[1] }

    versions        =   ch_versions
    

}