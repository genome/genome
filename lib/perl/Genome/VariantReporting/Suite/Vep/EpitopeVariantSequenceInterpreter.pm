package Genome::VariantReporting::Suite::Vep::EpitopeVariantSequenceInterpreter;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw(uniq each_array);

class Genome::VariantReporting::Suite::Vep::EpitopeVariantSequenceInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    doc => 'Output wildtype and mutant variant sequences for epitope prediction',
    has => {
        peptide_sequence_length => {
            is => 'Number',
            doc => 'The length of the peptide sequences',
            valid_values => [17, 21, 31],
            default_value => 21,
        },
    },
};

sub name {
    return 'epitope-variant-sequence';
}

sub requires_annotations {
    return ('vep');
}

sub field_descriptions {
    return (
        variant_sequences => "A hash of wildtype and mutant variant sequences for the entry's transcripts",
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my $vep_parser = new Genome::File::Vcf::VepConsequenceParser($entry->{header});

    my %return_values;
    ALLELE: for my $variant_allele (@$passed_alt_alleles) {
        my @transcripts = $vep_parser->transcripts($entry, $variant_allele);

        TRANSCRIPT: for my $transcript (@transcripts) {
            my @consequences = uniq map {split /\&/, lc($_)} grep {defined($_)} $transcript->{consequence};
            if (grep {$_ eq 'missense_variant'} @consequences) {
                my $position = $transcript->{protein_position} - 1;
                my ($wildtype_amino_acid, $mutant_amino_acid) = split('/', $transcript->{amino_acids});
                my @full_wildtype_sequence = split('', $transcript->{wildtypeprotein});
                if ($wildtype_amino_acid ne $full_wildtype_sequence[$position]) {
                    next TRANSCRIPT;
                }
                else {
                    my ($mutation_position, @wildtype_subsequence) = $self->get_wildtype_subsequence($position, \@full_wildtype_sequence, $entry);
                    my @mutant_subsequence = @wildtype_subsequence;
                    $mutant_subsequence[$mutation_position] = $mutant_amino_acid;
                    my @designations = qw(WT MT);
                    my @subsequences = (\@wildtype_subsequence, \@mutant_subsequence);
                    my $iterator = each_array( @designations, @subsequences );
                    while ( my ($designation, $subsequence) = $iterator->() ) {
                        my $header = ">$designation." . $transcript->{symbol} . '.p.' . $wildtype_amino_acid . $transcript->{protein_position} . $mutant_amino_acid;
                        my $sequence = join('', @$subsequence);
                        $return_values{$variant_allele}->{variant_sequences}->{$header} = $sequence;
                    }
                }
            }
            else {
                #Die?
            }
        }
    }

    return %return_values;
}
sub get_wildtype_subsequence {
    my ($self, $position, $full_wildtype_sequence_ref, $entry) = @_;
    my @full_wildtype_sequence = @$full_wildtype_sequence_ref;

    my $peptide_sequence_length = $self->peptide_sequence_length;
    #If the wildtype sequence is shorter than the desired peptide sequence
    #length we use the wildtype sequence length instead so that the extraction
    #algorithm below works correctly
    if (scalar(@full_wildtype_sequence) < $peptide_sequence_length) {
        $peptide_sequence_length = scalar(@full_wildtype_sequence);
        $self->status_message(
            'Wildtype sequence length is shorter than desired peptide sequence length at position (%s, %s). Using wildtype sequence length (%s) instead.',
            $entry->{chrom},
            $entry->{position},
            $peptide_sequence_length
        );
    }

    # We want to extract a subset from @full_wildtype_sequence that is
    # $peptide_sequence_length long so that the $position ends
    # up in the middle of the extracted sequence.
    # If the $position is too far toward the beginning or end of
    # @full_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.
    my (@wildtype_subsequence, $mutation_position);
    my $one_flanking_sequence_length = ($peptide_sequence_length - 1) / 2;
    if (distance_from_start($position, @full_wildtype_sequence) < $one_flanking_sequence_length) {
        @wildtype_subsequence = @full_wildtype_sequence[ 0 ... ($peptide_sequence_length - 1) ];
        $mutation_position = $position;
    }
    elsif (distance_from_end($position, @full_wildtype_sequence) < $one_flanking_sequence_length) {
        @wildtype_subsequence = @full_wildtype_sequence[ ($#full_wildtype_sequence - $peptide_sequence_length + 1) ... $#full_wildtype_sequence];
        $mutation_position = $peptide_sequence_length - distance_from_end($position, @full_wildtype_sequence) - 1;
    }
    elsif (
        (distance_from_start($position, @full_wildtype_sequence) >= $one_flanking_sequence_length) &&
        (distance_from_end($position, @full_wildtype_sequence) >= $one_flanking_sequence_length)
    ) {
        @wildtype_subsequence = @full_wildtype_sequence[ ($position - $one_flanking_sequence_length) ... ($position + $one_flanking_sequence_length) ];
        $mutation_position = ($peptide_sequence_length - 1) / 2;
    }
    else {
        die $self->error_message(
            'Something went wrong during the retrieval of the wildtype sequence at position (%s, %s, %s)',
            $entry->{chrom},
            $entry->{position},
            join('', @full_wildtype_sequence)
        );
    }

    return $mutation_position, @wildtype_subsequence;
}

#This subroutine is a bit funky but it was designed that way to mirror
#distance_from_end to increase code readability from the caller's perspective
sub distance_from_start {
    my $position = shift;
    my @array = @_;

    return $position;
}

sub distance_from_end {
    my $position = shift;
    my @array = @_;

    return $#array - $position;
}


1;
