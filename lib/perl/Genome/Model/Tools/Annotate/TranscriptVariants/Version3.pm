package Genome::Model::Tools::Annotate::TranscriptVariants::Version3;

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

=pod
=head1 Name

Genome::Transcript::VariantAnnotator

=head1 Synopsis

Given a variant, all transcripts affected by that variant are annotated and returned

=head1 Usage

# Variant file tab delimited, columns are chromosome, start, stop, reference, variant
# Need to infer variant type (SNP, DNP, INS, DEL) as well
my $variant_file = variants.tsv;
my @headers = qw/ chromosome_name start stop reference variant /;
my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $variant_file,
        headers => \@headers,
        separator => "\t",
        is_regex => 1,
        );

my $model = Genome::Model->get(name => 'NCBI-human.combined-annotation');
my $build = $model->build_by_version('54_36p');
my $iterator = $build->transcript_iterator;
my $window = Genome::Utility::Window::Transcript->create(
        iterator => $iterator,
        range => 50000,
        );
my $annotator = Genome::Transcript::VariantAnnotator->create(
        transcript_window => $window
        );

while (my $variant = $reader->next) {
    my @annotations = annotator->transcripts($variant);
}

=head1 Methods

=head2 transcripts

=over

=item I<Synopsis>   gets all annotations for a variant

=item I<Arguments>  variant (hash; see 'Variant Properites' below)

=item I<Returns>    annotations (array of hash refs; see 'Annotation' below)

=back

=head2 prioritized_transcripts

=over

=item I<Synopsis>   Gets one prioritized annotation per gene for a variant(snp or indel)

    =item I<Arguments>  variant (hash; see 'Variant properties' below)

    =item I<Returns>    annotations (array of hash refs; see 'Annotation' below)

    =back

    =head2 prioritized_transcript

    =over

    =item I<Snynopsis>  Gets the highest priority transcript affected by variant

    =item I<Arguments>  variant (hash, see 'Variant properties' below)

    =item I<Returns>    annotations (array of hash refs; see 'Annotation' below)

    =back

    =head1 Variant Properties

    =over

    =item I<chromosome_name>  The chromosome of the variant

    =item I<start>            The start position of the variant

    =item I<stop>             The stop position of the variant

    =item I<variant>          The snp base

    =item I<reference>        The reference base at the position

    =item I<type>             snp, dnp, ins, or del

    =back

    =head1 Annotation Properties

    =over

    =item I<transcript_name>    Name of the transcript

    =item I<transcript_source>  Source of the transcript

    =item I<strand>             Strand of the transcript

    =item I<c_position>         Relative position of the variant

    =item I<trv_type>           Called Classification of variant

=item I<priority>           Priority of the trv_type (only from get_prioritized_annotations)

    =item I<gene_name>          Gene name of the transcript

    =item I<intensity>          Gene intenstiy

    =item I<detection>          Gene detection

    =item I<amino_acid_length>  Amino acid length of the protein

    =item I<amino_acid_change>  Resultant change in amino acid in snp is in cds_exon

    =item I<variations>         Hashref w/ keys of known variations at the variant position

    =item I<type>               snp, ins, or del

    =back

    =head1 See Also

    B<Genome::Model::Command::Report>

    =head1 Disclaimer

    Copyright (C) 2008 Washington University Genome Sequencing Center

    This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

    Core Logic:

    B<Xiaoqi Shi> I<xshi@genome.wustl.edu>

    Optimization:

    B<Eddie Belter> I<ebelter@watson.wustl.edu>

    B<Gabe Sanderson> l<gsanders@genome.wustl.edu>

    B<Adam Dukes l<adukes@genome.wustl.edu>

    B<Brian Derickson l<bdericks@genome.wustl.edu>

    =cut

#$HeadURL$
#$Id$
