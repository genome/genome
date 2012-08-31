package Genome::Model::Tools::Annotate::AminoAcidChange;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::AminoAcidChange{
    is => 'Genome::Model::Tools::Annotate',
    has => [
        amino_acid_change_string => {
            is => 'Text',
            is_input => 1,
            doc => 'Amino Acid change string to check', 
        },
    ],
    has_optional => [
        amino_acid_change_results => {
            is => 'SCALAR', 
            is_output => 1,
            doc => "This is populated with an array reference containing the results of the amino acid change check post execute",
            is_optional => 1,
        },
    ],
};


sub help_synopsis {
    return <<EOS
my \$lookup = Genome::Model::Tools::Annotate::AminoAcidChange->execute(amino_acid_change_string => \$amino_acid_change_string);
my \$results_array_ref = \$lookup->amino_acid_change_results;
EOS
}

sub execute {
    my $self = shift;

    my $results = 
        eval { check_amino_acid_change_string(
                   amino_acid_change_string => $self->amino_acid_change_string)
        };
    if ($@) {
        $self->error_message($@);
        return 0;
    } else {
        $self->amino_acid_change_results($results);
        return 1;
    }
}

sub check_amino_acid_change_string{
    my %args = @_;
    my $amino_acid_change_string = $args{'amino_acid_change_string'};
    my ($residue1, $res_start, $residue2, $res_stop, $new_residue);

    unless (defined($amino_acid_change_string)) {
        return [$residue1, $res_start, $residue2, $res_stop, $new_residue];
    }

    #__VALIDATE
    $amino_acid_change_string =~ s/^p\.//x;
    if ($amino_acid_change_string =~ /^ (\D+) (\d+) _ (\D+) (\d+) (.*) $/x ) {
        ($residue1, $res_start, $residue2, $res_stop, $new_residue) =
            ($1, $2, $3, $4, $5);
    } elsif ($amino_acid_change_string =~ /^ (\D+) (\d+) (\D+) (.*) $/x ) {
        ($residue1, $res_start, $residue2, $new_residue) =
            ($1, $2, $3, $4);
        $res_stop = $res_start;
    } elsif ($amino_acid_change_string =~ /^ (\d+) (.*) $/x ) {
        ($res_start, $residue2, $res_stop, $new_residue) =
            ($1, $2);
        $residue1 = '*';
        $res_stop = $res_start;
        $new_residue = $residue2;
    }
    if (defined($new_residue)) {
        $new_residue =~ s/^ > //x;
    }
    $new_residue ||= '';

    return [$residue1, $res_start, $residue2, $res_stop, $new_residue];
}

1;
