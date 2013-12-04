package Genome::Site::TGI::Synchronize::Classes::LibrarySummary; 

use strict;
use warnings;

=pod
FULL_NAME            VARCHAR2 (64)                    {null} NOT NULL ok [name]
GENOTYPE_QC_SEQ_ID   NUMBER   (10)                    {null} {null}   NOT_SYNCED
LIBRARY_ID           NUMBER   (20)                    {null} NOT NULL ok [id]
LIBRARY_INSERT_SIZE  VARCHAR2 (64)                    {null} {null}   ok
ORIGINAL_INSERT_SIZE VARCHAR2 (64)                    {null} {null}   ok
PROTOCOL             VARCHAR2 (64)                    {null} {null}   ok
TRANSCRIPT_STRAND    VARCHAR2 (16)                    {null} {null}   ok
SAMPLE_ID            NUMBER   (20)                    {null} NOT NULL ok
=cut

class Genome::Site::TGI::Synchronize::Classes::LibrarySummary {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'LIBRARY_SUMMARY',
    id_by => [
        id                      => { is => 'Number', column_name => 'LIBRARY_ID', },
    ],
    has => [
        name                    => { is => 'Text', column_name => 'FULL_NAME', },
        sample_id               => { is => 'Number', },
    ],
    has_optional => [
        protocol                => { is => 'Text', },
        library_insert_size     => { is => 'Text', },
        original_insert_size    => { is => 'Text', },
        transcript_strand       => { is => 'Text', },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub entity_name { return 'library' }

sub properties_to_copy {# 6
    return ( 'id', 'name', 'sample_id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {# 3
    return (qw/ 
        protocol
        library_insert_size
        original_insert_size
        transcript_strand
        /);
}

sub lims_property_name_to_genome_property_name {
    my ($class, $name) = @_;
    my %lims_to_genome = (
        full_name => 'name',
    );
    return $lims_to_genome{$name} if exists $lims_to_genome{$name};
    return $name;
}



1;

