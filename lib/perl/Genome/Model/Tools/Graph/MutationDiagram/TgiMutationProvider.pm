package Genome::Model::Tools::Graph::MutationDiagram::TgiMutationProvider;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Graph::MutationDiagram::TgiMutationProvider {
    is => 'Genome::Model::Tools::Graph::MutationDiagram::MutationProvider',
    has => [
        input_file => {
            is => 'File',
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
        headers => [qw(na0 na1 na2 na3 na4 na5 hugo transcript_name na8 na9 na10 na11 na12 class c_position aa_change)],
        ignore_extra_columns => 1,
    );
    $self->reader($reader);
    return $self;
}

sub next {
    my $self = shift;
    while (my $data = $self->reader->next) {
        my %params;
        my $transcript_name = $data->{transcript_name};
        next unless defined($transcript_name) && $transcript_name !~ /^\s*$/ && $transcript_name ne '-';
        my ($residue1, $res_start, $residue2, $res_stop, $new_residue) = @{Genome::Model::Tools::Annotate::AminoAcidChange::check_amino_acid_change_string(amino_acid_change_string => $data->{aa_change})};
        my $mutation = $data->{aa_change};
        $mutation =~ s/p\.//g;

        if($data->{aa_change} =~ /^e/ && $data->{c_position} !~ /NULL/) {
#this is a splice site mutation
            $res_start = c_position_to_amino_acid_pos($data->{c_position});
        }
        $params{mutation} = $mutation;
        $params{hugo} = $data->{hugo};
        $params{transcript_name} = $transcript_name;
        $params{protein_position} = $res_start;
        $params{class} = $data->{class};
        return \%params;
    }
    return;
}

sub c_position_to_amino_acid_pos {
    my ($c_position) = @_;
    ($c_position) = split(/[\+\-]/, $c_position); #ignore +1 etc and change info
    $c_position =~ s/[^0-9]//g; #strip off c.
    my $aa_pos = sprintf("%d", ($c_position / 3));
    return $aa_pos;
}

1;

