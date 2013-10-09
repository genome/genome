package Genome::Model::Tools::Annotate::TranscriptVariants::Version4;

use strict;
use warnings;

use Data::Dumper;
use Genome;
use File::Temp;
use List::Util qw/ max min /;
use List::MoreUtils qw/ uniq /;
use Bio::Seq;
use Bio::Tools::CodonTable;
use DateTime;
use Carp;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Genome::Model::Tools::Annotate::TranscriptVariants::Base',

    has => [
        reference_sequence_id => {
            is => "Text",
        },
        eids => {
            is_transient => 1,
            is_optional => 1,
            doc => "Temporary variable used for intermediate calculation.",
        },
    ],

    doc => q(Do proper intersections between variations and transcript structures by considering both entities' start and stop positions rather than just the start position of the variation.),

);


sub transcript_status_priorities {
    return (
        known               => 1,
        putative            => 2,
        novel               => 3,
    );
}


sub is_mitochondrial {
    my ($self, $chrom_name) = @_;

    #we use the mitochondrial_codon_translator if the chromosome is either the M or MT.  Everything else should use the normal translator
    return $chrom_name =~ /^MT?/;
}


sub cache_gene_names {
    my $self = shift;

    if (!defined $self->{_cached_chromosome}) {
        Genome::ExternalGeneId->get(
            data_directory => $self->data_directory,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl",
        );
        Genome::ExternalGeneId->get(
            data_directory => $self->data_directory,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl_default_external_name",
        );
        Genome::ExternalGeneId->get(
            data_directory => $self->data_directory,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl_default_external_name_db",
        );
    }
}

sub filter_and_partition_structures {
}

sub specialized_deletion_annotation {
}

# Annotates a single transcript-substructure/variant pair
sub _transcript_substruct_annotation {
    my ($self, $substruct, %variant) = @_;

    my %result = $self->SUPER::_transcript_substruct_annotation($substruct, %variant);

    my $dumper_string = $substruct->id;
    my ($default_gene_name, $ensembl_gene_id, $gene_name_source);
    unless ($self->eids) {
        my %new;
        $self->eids(\%new);
    }
    if ($self->eids and $self->eids->{$substruct->transcript_gene_id}) {
        ($default_gene_name,$ensembl_gene_id,$gene_name_source) = split(/,/, $self->eids->{$substruct->transcript_gene_id});
    }
    else {
        my @e1 = Genome::ExternalGeneId->get(data_directory => $substruct->data_directory,
            gene_id => $substruct->transcript_gene_id,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl_default_external_name");
        if ($e1[0]) {
            $default_gene_name = $e1[0]->id_value;
        }
        unless ($default_gene_name) {
            $self->warning_message("Ensembl gene name missing for substruct: $dumper_string");
            $default_gene_name = "Unknown";
        }
        my @e2 = Genome::ExternalGeneId->get(data_directory => $substruct->data_directory,
            gene_id => $substruct->transcript_gene_id,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl");
        if ($e2[0]) {
            $ensembl_gene_id = $e2[0]->id_value;
        }
        unless ($ensembl_gene_id) {
            $self->warning_message("Ensembl stable gene id missing for substruct: $dumper_string");
            $ensembl_gene_id = "Unknown";
        }
        my @e3 = Genome::ExternalGeneId->get(data_directory => $substruct->data_directory,
            gene_id => $substruct->transcript_gene_id,
            reference_build_id => $self->reference_sequence_id,
            id_type => "ensembl_default_external_name_db");
        if ($e3[0]) {
            $gene_name_source = $e3[0]->id_value;
        }
        unless ($gene_name_source) {
            $self->warning_message("Ensembl gene name source missing for substruct: $dumper_string");
            $gene_name_source = "Unknown";
        }
        $self->eids->{$substruct->transcript_gene_id} = join(',',$default_gene_name, $ensembl_gene_id, $gene_name_source);
    }

    $result{default_gene_name} = $default_gene_name;
    $result{gene_name_source} = $gene_name_source;
    $result{ensembl_gene_id} = $ensembl_gene_id;

    return %result;
}

sub should_update_variant_attributes {
    my ($self, $variant, $structure) = @_;

    unless ($self->{'get_frame_shift_sequence'}) {
        # If we're inspecting the entire sequence, don't chop the variant down...
        if ($variant->{stop} > $structure->structure_stop and $variant->{type} eq 'DEL') {
            return 1;
        }
    }

    return;
}

sub get_dnp_snp_trv_type {
    my ($self, $original_aa, $mutated_aa) = @_;

    if (!defined $mutated_aa or !defined $original_aa) {
        return 'silent';
    }
    elsif ($mutated_aa eq $original_aa) {
        return 'silent';
    }
    else {
        my ($reduced_original_aa, $reduced_mutated_aa, $offset) = $self->_reduce(
            $original_aa, $mutated_aa);

        if (index($reduced_mutated_aa, '*') != -1) {
            return 'nonsense';
        }
        elsif (index($reduced_original_aa, '*') != -1) {
            return 'nonstop';
        }
        else {
            return 'missense';
        }
    }
}

sub get_dnp_snp_protein_data {
    my ($self, $original_aa, $mutated_aa, $protein_position) = @_;

    if (!defined $mutated_aa or !defined $original_aa) {
        return ("NULL", $protein_position);
    }
    elsif ($mutated_aa eq $original_aa) {
        return ("p." . $original_aa . $protein_position, $protein_position);
    }
    else {
        my ($reduced_original_aa, $reduced_mutated_aa, $offset) = $self->_reduce(
            $original_aa, $mutated_aa);
        $protein_position += $offset;

        return ("p." . $reduced_original_aa . $protein_position . $reduced_mutated_aa,
            $protein_position);
    }
}

# Taken from Genome::Transcript
# Given a version and species, find the imported reference sequence build
sub get_reference_build_for_transcript {
    my($self, $structure) = @_;

    my ($version) = $structure->transcript_version =~ /^\d+_(\d+)[a-z]/;
    my $species = $structure->transcript_species;

    unless ($self->{'_reference_builds'}->{$version}->{$species}) {
        my $build = Genome::Model::Build->get($self->reference_sequence_id);
        confess "Could not get build version $version" unless $build;

        $self->{'_reference_build'}->{$version}->{$species} = $build;
    }
    return $self->{_reference_build}->{$version}->{$species};
}


sub bound_relative_stop {
    my ($self, $relative_stop, $limit) = @_;

    #it is possible that the variant goes off the end of the transcript.  In this case,
    #we need to adjust the relative stop.

    return min($relative_stop, $limit);
}

1;
