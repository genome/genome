package Genome::Model::Tools::Graph::MutationDiagram::VepMutationProvider;

use strict;
use warnings;
use Genome;

my @VEP_MUTATION_PRIORITY = (
    'ESSENTIAL_SPLICE_SITE',
    'FRAMESHIFT_CODING',
    'STOP_GAINED',
    'NON_SYNONYMOUS_CODING'
);

my %VEP_MUTATION_PRIORITIES;
@VEP_MUTATION_PRIORITIES{@VEP_MUTATION_PRIORITY} = 0..$#VEP_MUTATION_PRIORITY;

class Genome::Model::Tools::Graph::MutationDiagram::VepMutationProvider {
    is => 'Genome::Model::Tools::Graph::MutationDiagram::MutationProvider',
    has => [
        input_file => {
            is => 'File',
        },
        vep_frequency_field => {
            is => 'Text',
        },
    ],
    has_transient_optional => [
        reader => {
            is => 'Genome::Utility::IO::SeparatedValueReader',
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->input_file,
        separator => "\t",
        headers => [qw(something0 something1 something2 hugo transcript_name something5 class something7 something8 protein_pos aa_change extra)],
        ignore_extra_columns => 1,
    );
    my $header = $reader->next;
    $self->reader($reader);
    return $self;
}

sub next {
    my $self = shift;
    while (my $data = $self->reader->next) {
        my %params;
        next unless(defined($data->{transcript_name}) && $data->{transcript_name} !~ /^\s*$/);
        $params{hugo} = $data->{hugo};
        $params{transcript_name} = $data->{transcript_name};
        $params{protein_position} = $data->{protein_pos};
        $params{mutation} = $self->_calculate_mutation($data->{aa_change}, $data->{protein_pos});
        next if $params{mutation} eq '--';
        $params{class} = $self->_get_vep_mutation_class($data->{class});
        $params{frequency} = $self->_calculate_frequency($data->{extra});
        return \%params;
    }
    return;
}

sub _calculate_mutation {
    my $self = shift;
    my $aa_change = shift;
    my $protein_pos = shift;
    my ($orig_aa, $new_aa) = split("/", $aa_change);
    $orig_aa |= '';
    $new_aa |= '';
    my $mutation = join($protein_pos, $orig_aa, $new_aa);
    return $mutation;
}

sub _get_vep_mutation_class {
    my $self = shift;
    my $type = shift;
    my @types = split(",", $type);
    my @type_matches = grep {defined $VEP_MUTATION_PRIORITIES{$_}} @types;
    return $type unless @type_matches;
    my @priorities = @VEP_MUTATION_PRIORITIES{@type_matches};
    my $idx = argmin(@priorities);
    return $type_matches[$idx];
}

sub _calculate_frequency {
    my $self = shift;
    my $extra_fields = shift;
    my $extra = $self->_get_vep_extra_fields_hash($extra_fields);
    my $frequency = 1;
    $frequency = $extra->{$self->vep_frequency_field} if exists $extra->{$self->vep_frequency_field};
    return $frequency;
}

sub _get_vep_extra_fields_hash {
    my $self = shift;
    my $extra = shift;
    return { map { split("=") } split(";", $extra) }
}

1;

