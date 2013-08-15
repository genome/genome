package Genome::Site::TGI::Synchronize::Classes::RegionIndex454;

use strict;
use warnings;

use Genome;

=comment
# RUN_REGION_454
REGION_ID                NUMBER   (12)                    {null} NOT NULL NOT_SYNCED [from ri454]
ANALYSIS_NAME            VARCHAR2 (255)                   {null} NOT NULL NOT_SYNCED
REGION_NUMBER            NUMBER   (12)                    {null} NOT NULL ok
TOTAL_RAW_WELLS          NUMBER   (12)                    {null} NOT NULL NOT_SYNCED
TOTAL_KEY_PASS           NUMBER   (12)                    {null} NOT NULL NOT_SYNCED
INCOMING_DNA_NAME        VARCHAR2 (64)                    {null} NOT NULL NOT_SYNCED
COPIES_PER_BEAD          NUMBER   (16,4)                  {null} NOT NULL NOT_SYNCED
RUN_NAME                 VARCHAR2 (255)                   {null} NOT NULL ok
KEY_PASS_WELLS           NUMBER   (12)                    {null} NOT NULL NOT_SYNCED
PREDICTED_RECOVERY_BEADS NUMBER   (16,4)                  {null} {null}   NOT_SYNCED
FC_ID                    NUMBER   (10)                    {null} {null}   NOT_SYNCED
SAMPLE_SET               VARCHAR2 (64)                    {null} {null}   NOT_SYNCED
RESEARCH_PROJECT         VARCHAR2 (64)                    {null} {null}   NOT_SYNCED
PAIRED_END               NUMBER   (1)                     {null} {null}   ok
SAMPLE_NAME              VARCHAR2 (64)                    {null} {null}   NOT_SYNCED
LIBRARY_NAME             VARCHAR2 (64)                    {null} {null}   NOT_SYNCED
BEADS_LOADED             NUMBER   (12)                    {null} {null}   NOT_SYNCED
SS_ID                    NUMBER   (15)                    {null} {null}   NOT_SYNCED
SUPERNATANT_BEADS        NUMBER   (8)                     {null} {null}   NOT_SYNCED
SAMPLE_ID                NUMBER   (20)                    {null} {null}   NOT_SYNCED
LIBRARY_ID               NUMBER   (20)                    {null} {null}   NOT_SYNCED [from ri454]
# REGION_INDEX_454
SEQ_ID             NUMBER   (15)                    {null} NOT NULL ok [id]
REGION_ID          NUMBER   (12)                    {null} NOT NULL ok
INDEX_SEQUENCE     VARCHAR2 (255)                   {null} {null}   ok
NUM_READS          NUMBER   (12)                    {null} NOT NULL ok [total_reads]
NUM_BASES          NUMBER   (12)                    {null} NOT NULL ok [total_bases_read]
TAGGED_DNA_ID      NUMBER   (20)                    {null} {null}   NOT_SYNCED
CREATION_EVENT_ID  NUMBER   (15)                    {null} {null}   NOT_SYNCED
LIBRARY_ID         NUMBER   (20)                    {null} {null}   ok
LIBRARY_STRATEGY   VARCHAR2 (30)                    {null} {null}   NOT_SYNCED
EXPERIMENT_PURPOSE VARCHAR2 (32)                    {null} {null}   NOT_SYNCED
=cut

class Genome::Site::TGI::Synchronize::Classes::RegionIndex454 {
    table_name => <<'SQL'
        (
            select 
                --ri454
                to_char(ri454.seq_id) id,
                ri454.index_sequence index_sequence,
                ri454.library_id,
                ri454.region_id,
                ri454.num_bases total_bases_read,
                ri454.num_reads total_reads,
                --rr454
                rr454.paired_end is_paired_end,
                rr454.region_number,
                rr454.run_name run_name,
                ( 
                 case when ri454.index_sequence is null
                  then to_char(rr454.region_number)
                  else to_char(rr454.region_number) || '-' || ri454.index_sequence 
                 end
                ) subset_name,
                --constants
                '454' sequencing_platform
            from region_index_454 ri454
            join run_region_454 rr454 on rr454.region_id = ri454.region_id
        ) region_index_454
SQL
    ,
    id_by => [
        id => { is => 'Text', },
    ],
    has_optional => [
        index_sequence => { is => 'Text', },
        is_paired_end => { is => 'Text', },
        library_id => { is =>'Text', },
        region_id => { is => 'Text', },
        region_number => { is => 'Text', },
        run_name => { is => 'Text', },
        sequencing_platform => { is => 'Text', },
        subset_name => { is => 'Text', },
        total_reads => { is => 'Text', },
        total_bases_read => { is => 'Text', },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub properties_to_copy {# 9
    return ( 'id', 'library_id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {# 7
    return (qw/ 
        index_sequence
        is_paired_end
        region_id
        region_number
        run_name
        total_reads
        total_bases_read
        /);
}

sub lims_property_name_to_genome_property_name {
    my ($class, $name) = @_;
    my %lims_to_genome = (
        num_reads => 'total_reads',
        num_bases => 'total_bases_read',
        paired_end => 'is_paired_end',
    );
    return $lims_to_genome{$name} if exists $lims_to_genome{$name};
    return $name;
}


