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
    my %transcript_count;
    ALLELE: for my $variant_allele (@$passed_alt_alleles) {
        my @transcripts = $vep_parser->transcripts($entry, $variant_allele);

        TRANSCRIPT: for my $transcript (@transcripts) {
            $transcript_count{$transcript->{feature}}++;
            my @consequences = uniq map {split /\&/, lc($_)} grep {defined($_)} $transcript->{consequence};
            my ($wildtype_amino_acid, $mutant_amino_acid) = split('/', $transcript->{amino_acids});
            my $full_wildtype_sequence = $transcript->{wildtypeprotein};
            my ($position, $wildtype_amino_acid_length);
            if (grep {$_ eq 'missense_variant' || $_ eq 'inframe_insertion'} @consequences) {
                if ($wildtype_amino_acid eq '-') {
                    ($position) = split('-', $transcript->{protein_position});
                    $wildtype_amino_acid_length = 0;
                }
                else {
                    $position = $transcript->{protein_position} - 1;
                    $wildtype_amino_acid_length = length($wildtype_amino_acid);
                }
            }
            elsif (grep {$_ eq 'inframe_deletion'} @consequences) {
                if ($mutant_amino_acid eq '-') {
                    $position = $transcript->{protein_position} - 1;
                    $wildtype_amino_acid_length = length($wildtype_amino_acid);
                    $mutant_amino_acid = '';
                }
                else {
                    $position = (split('-', $transcript->{protein_position}))[0] - 1;
                    $wildtype_amino_acid_length = length($wildtype_amino_acid);
                }
            }
            elsif (grep {$_ eq 'frameshift_variant'} @consequences) {
                $position = $transcript->{protein_position} - 1;
            }
            else {
                next TRANSCRIPT;
            }

            next TRANSCRIPT if position_out_of_bounds($position, $full_wildtype_sequence);

            my ($wildtype_subsequence, $mutant_subsequence, $variant_id);
            if (grep {$_ eq 'frameshift_variant'} @consequences) {
                ($wildtype_subsequence, $mutant_subsequence) = $self->get_frameshift_subsequences($position, $full_wildtype_sequence, $entry);
                $mutant_subsequence .= $transcript->{downstreamprotein};
                $variant_id = $transcript->{symbol} . '_' . $transcript->{feature} . '_' . $transcript_count{$transcript->{feature}} . '.FS.' . $transcript->{protein_position};
            }
            else {
                my $mutation_position;
                ($mutation_position, $wildtype_subsequence) = $self->get_wildtype_subsequence($position, $full_wildtype_sequence, $entry, $wildtype_amino_acid_length);
                $mutant_subsequence = $wildtype_subsequence;
                substr($mutant_subsequence, $mutation_position, $wildtype_amino_acid_length) = $mutant_amino_acid;
                $variant_id = $transcript->{symbol} . '_' . $transcript->{feature} . '_' . $transcript_count{$transcript->{feature}} . '.' . $wildtype_amino_acid . $transcript->{protein_position} . ($mutant_amino_acid || '-');
            }
            my @designations = qw(WT MT);
            my @subsequences = ($wildtype_subsequence, $mutant_subsequence);
            my $iterator = each_array( @designations, @subsequences );
            while ( my ($designation, $subsequence) = $iterator->() ) {
                my $header = ">$designation." . $variant_id;
                $return_values{$variant_allele}->{variant_sequences}->{$header} = $subsequence;
            }
        }

        if (!defined($return_values{$variant_allele})) {
            $return_values{$variant_allele}->{variant_sequences} = '';
        }
    }

    return %return_values;
}

sub determine_peptide_sequence_length {
    my ($self, $full_wildtype_sequence_length, $entry) = @_;

    my $peptide_sequence_length = $self->peptide_sequence_length;

    #If the wildtype sequence is shorter than the desired peptide sequence
    #length we use the wildtype sequence length instead so that the extraction
    #algorithm below works correctly
    if ($full_wildtype_sequence_length < $peptide_sequence_length) {
        $peptide_sequence_length = $full_wildtype_sequence_length;
        $self->status_message(
            'Wildtype sequence length is shorter than desired peptide sequence length at position (%s, %s). Using wildtype sequence length (%s) instead.',
            $entry->{chrom},
            $entry->{position},
            $peptide_sequence_length
        );
    }

    return $peptide_sequence_length;
}

sub determine_flanking_sequence_length {
    my ($self, $full_wildtype_sequence_length, $entry) = @_;
    return ($self->determine_peptide_sequence_length($full_wildtype_sequence_length, $entry) - 1) / 2;
}

sub get_frameshift_subsequences {
    my ($self, $position, $full_wildtype_sequence, $entry) = @_;

    my $one_flanking_sequence_length = $self->determine_flanking_sequence_length(length($full_wildtype_sequence), $entry);
    my $start_position = $position - $one_flanking_sequence_length;
    my $wildtype_subsequence_length =  2 * $one_flanking_sequence_length;
    my $mutation_subsequence_length =  $one_flanking_sequence_length;
    if ($start_position < 0) {
        $start_position = 0;
        $wildtype_subsequence_length = $position + $one_flanking_sequence_length;
        $mutation_subsequence_length = $position;
    }
    my $wildtype_subsequence = substr($full_wildtype_sequence, $start_position, $wildtype_subsequence_length);
    my $mutation_start_subsequence = substr($full_wildtype_sequence, $start_position, $mutation_subsequence_length);
    return ($wildtype_subsequence, $mutation_start_subsequence);
}

sub get_wildtype_subsequence {
    my ($self, $position, $full_wildtype_sequence, $entry, $wildtype_amino_acid_length) = @_;

    my $one_flanking_sequence_length = $self->determine_flanking_sequence_length(length($full_wildtype_sequence), $entry);
    my $peptide_sequence_length = 2 * $one_flanking_sequence_length + $wildtype_amino_acid_length;

    # We want to extract a subset from $full_wildtype_sequence that is
    # $peptide_sequence_length long so that the $position ends
    # up in the middle of the extracted sequence.
    # If the $position is too far toward the beginning or end of
    # $full_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.
    my ($wildtype_subsequence, $mutation_position);
    if (distance_from_start($position, $full_wildtype_sequence) < $one_flanking_sequence_length) {
        $wildtype_subsequence = substr($full_wildtype_sequence, 0, $peptide_sequence_length);
        $mutation_position = $position;
    }
    elsif (distance_from_end($position, $full_wildtype_sequence) < $one_flanking_sequence_length) {
        $wildtype_subsequence = substr($full_wildtype_sequence, (length($full_wildtype_sequence) - $peptide_sequence_length), $peptide_sequence_length);
        $mutation_position = $peptide_sequence_length - distance_from_end($position, $full_wildtype_sequence) - 1;
    }
    elsif (
        (distance_from_start($position, $full_wildtype_sequence) >= $one_flanking_sequence_length) &&
        (distance_from_end($position, $full_wildtype_sequence) >= $one_flanking_sequence_length)
    ) {
        $wildtype_subsequence = substr($full_wildtype_sequence, ($position - $one_flanking_sequence_length), $peptide_sequence_length);
        $mutation_position = $one_flanking_sequence_length;
    }
    else {
        die $self->error_message(
            'Something went wrong during the retrieval of the wildtype sequence at position (%s, %s, %s)',
            $entry->{chrom},
            $entry->{position},
            $full_wildtype_sequence
        );
    }

    return $mutation_position, $wildtype_subsequence;
}

#This subroutine is a bit funky but it was designed that way to mirror
#distance_from_end to increase code readability from the caller's perspective
sub distance_from_start {
    my $position = shift;
    my $sequence = shift;

    return $position;
}

sub distance_from_end {
    my $position = shift;
    my $sequence = shift;

    return length($sequence) - $position - 1;
}

sub position_out_of_bounds {
    my $position = shift;
    my $sequence = shift;

    return ($position > length($sequence) - 1);
}


1;
