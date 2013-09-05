package Genome::Model::Tools::Annotate::TranscriptVariants::Version2;

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

    doc => q(Do proper intersections between variations and transcript structures by considering both entities' start and stop positions rather than just the start position of the variation),

);


sub transcript_status_priorities {
    return (
        reviewed            => 1,
        validated           => 2,
        provisional         => 3,
        predicted           => 4,
        putative            => 4,
        model               => 5,
        inferred            => 6,
        known               => 7,
        annotated           => 8,
        known_by_projection => 9,
        novel               => 10,
        unknown             => 11,
    );
}


sub is_mitochondrial {
    my ($self, $chrom_name) = @_;

    #we use the mitochondrial_codon_translator if the chromosome is either the M or MT.  Everything else should use the normal translator
    return $chrom_name =~ /^MT?/;
}


sub cache_gene_names {
    # Nothing to cache for Version 2
}

sub filter_and_partition_structures {
}

sub specialized_deletion_annotation {
}

sub reference_sequence_id {
    my $self = shift;

     return $self->build->reference_sequence_id;
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

        my $model = Genome::Model::ImportedReferenceSequence->get(name => "NCBI-$species");
        confess "Could not get imported reference sequence model for $species!" unless $model;

        my $build = $model->build_by_version($version . '-lite');
        unless ($build) {
            warn "no $version-lite reference, trying the full reference...";
        }
        confess "Could not get build version $version from $species imported reference sequence model!" unless $build;

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
