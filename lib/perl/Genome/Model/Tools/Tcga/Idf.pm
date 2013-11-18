package Genome::Model::Tools::Tcga::Idf;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Tcga::Idf {
    has => [
    ],
};

sub resolve_maf_protocol {
    my $self = shift;
    my $protocol_db = shift;

    unless (defined $protocol_db->{"mutation filtering annotation and curation"}) {
        $protocol_db->{"mutation filtering annotation and curation"} = [{name => "genome.wustl.edu:maf_creation:data_consolidation:01",
                                                                description => "Automatic and manual filtering and curation of variants"}];
    }
    return $protocol_db->{"mutation filtering annotation and curation"}->[0]->{name};
}

sub resolve_mapping_protocol {
    my $self = shift;
    my $processing_profile = shift;
    my $protocol_db = shift;

    my $name = "genome.wustl.edu:alignment:".$processing_profile->id.":01";
    my $description = $processing_profile->name;
    if (defined $protocol_db->{"sequence alignment"}){
        my $found = 0;
        for my $protocol (@{$protocol_db->{"variant calling"}}) {
            if ($protocol->{name} eq $name) {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            push @{$protocol_db->{"variant calling"}}, {name => $name, description => $description};
        }      
    }
    else {
        $protocol_db->{"sequence alignment"} = [{name => "genome.wustl.edu:DNA_sequencing:Illumina:01",
                                                                description => "Illumina sequencing by synthesis"}];
    }
    return $name;
}

sub resolve_library_protocol {
    my $self = shift;
    my $protocol_db = shift;

    unless (defined $protocol_db->{"library preparation"}){
        $protocol_db->{"library preparation"} = [{name => "genome.wustl.edu:DNA_extraction:Illumina_DNASeq:01",
                                                                description => "Illumina library prep"}];
    }
    return $protocol_db->{"library preparation"}->[0]->{name};
}

sub resolve_variants_protocol {
    my $self = shift;
    my $processing_profile = shift;
    my $protocol_db = shift;

    my $name = "genome.wustl.edu:variant_calling:".$processing_profile->id.":01";
    my $description = $processing_profile->name;
    if (defined $protocol_db->{"variant calling"}){
        my $found = 0;
        for my $protocol (@{$protocol_db->{"variant calling"}}) {
            if ($protocol->{name} eq $name) {
                $found = 1;
                last;
            }
        }
        unless ($found) {
            push @{$protocol_db->{"variant calling"}}, {name => $name, description => $description};
        }
    }
    else {
        $protocol_db->{"variant calling"} = [{name => $name, description => $description}];
    }
    return $name;
}

1;

