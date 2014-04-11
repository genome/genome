package Genome::Model::Tools::EpitopePrediction::GenerateVariantSequences;

use strict;
use warnings;

use Genome;
use Workflow;

my $DEFAULT_TRV_TYPE = 'missense';


class Genome::Model::Tools::EpitopePrediction::GenerateVariantSequences {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input file is a tab-separated (TSV) output file from gmt annotate variant-protein. For more info, gmt annotate variant-protein --help',
        },
        output_directory => {
            is => 'Text',
            doc => 'The output location',
        },
        trv_type => {
            is => 'Text',
            is_optional => 1,
            doc => 'The type of mutation you want to output eg missense,nonsense',
            # the current code only works on missense. Will need furthur development for other trv_types.
            valid_values => [$DEFAULT_TRV_TYPE],
            default_value => $DEFAULT_TRV_TYPE,
        },
        peptide_sequence_length => {
            is => 'Number',
            doc => 'The length of the peptide sequences',
            valid_values => [17, 21, 31],
            default_value => 21,
        },
    ],
    has_output => [
        output_file => {
            is => 'Text',
            doc => 'The output FASTA file to write 21mer sequences for wildtype(WT) and mutant(MT) proteins',
            calculate_from => ['output_directory'],
            calculate => q| return File::Spec->join($output_directory, "variant_sequences_unfiltered.fa");|,
        },
    ],
};


sub help_brief {
    "FOR NETMHC : Outputs a FASTA file for  for wildtype(WT) and mutant(MT) proteins 21-mer sequences for MHC Class I epitope prediction",
}



sub execute {
    my $self = shift;

    my $input_fh = Genome::Sys->open_file_for_reading($self->input_file);
    Genome::Sys->validate_file_for_writing($self->output_file);

    while (my $line = $input_fh->getline) {
        chomp $line;
        $line =~ s/[*]$//g;
        my @protein_arr =  split(/\t/, $line);
        next unless $protein_arr[13] eq $self->trv_type;
        if ( $protein_arr[15] =~ /^p.([A-Z])(\d+)([A-Z])/ ) {
            my ($wildtype_amino_acid, $position, $mutant_amino_acid) = ($1, $2 - 1, $3);
            my $wildtype_sequence = $self->get_wildtype_sequence(@protein_arr);
            my @arr_wildtype_sequence = split('',$wildtype_sequence);

            if ($wildtype_amino_acid ne $arr_wildtype_sequence[$position]) {
                next;
            }
            else {
                my @mutant_arr = my @wildtype_arr = $self->get_wildtype_subsequence_for_printing($position, @arr_wildtype_sequence);
                my $midpoint = ($self->peptide_sequence_length - 1) / 2;
                $mutant_arr[$midpoint] = $mutant_amino_acid;
                $self->print_to_output(\@wildtype_arr, \@mutant_arr, \@protein_arr, $position);
            }
        }
    }

    close($input_fh);

    return 1;
}

sub print_to_output {
    my $self = shift;
    my $wildtype_arr_ref = shift;
    my $mutant_arr_ref = shift;
    my $protein_arr_ref = shift;
    my $position = shift;

    if ($wildtype_arr_ref) {
        for my $designation ('WT', 'MT') {
            my $arr = $designation eq 'WT' ? $wildtype_arr_ref : $mutant_arr_ref;
            $self->print_fasta_entry(
                designation => $designation,
                arr         => $arr,
                protein_arr => $protein_arr_ref,
            );
        }
    }
    else {
        my $output_fh = Genome::Sys->open_file_for_appending($self->output_file);
        print $output_fh "NULL"."\t".$position."\n";
        close($output_fh);
    }
}

sub get_wildtype_subsequence_for_printing {
    my $self = shift;
    my $position = shift;
    my @arr_wildtype_sequence = @_;

    # We want to extract a subset from @arr_wildtype_sequence that is
    # $self->peptide_sequence_length long so that the $position ends
    # up in the middle of the subsequence.
    # If the $position is too far toward the beginning or end of
    # @arr_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.

    my @wildtype_arr;
    my $one_flanking_sequence_length = ($self->peptide_sequence_length - 1) / 2;
    if (distance_from_start($position, @arr_wildtype_sequence) < $one_flanking_sequence_length) {
        @wildtype_arr = @arr_wildtype_sequence[ 0 ... ($self->peptide_sequence_length - 1) ];
    }
    elsif (distance_from_end($position, @arr_wildtype_sequence) < $one_flanking_sequence_length) {
        @wildtype_arr = @arr_wildtype_sequence[ ($#arr_wildtype_sequence - $self->peptide_sequence_length) ... $#arr_wildtype_sequence];
    }
    elsif (
        (distance_from_start($position, @arr_wildtype_sequence) >= $one_flanking_sequence_length) &&
        (distance_from_end($position, @arr_wildtype_sequence) >= $one_flanking_sequence_length)
    ) {
        @wildtype_arr = @arr_wildtype_sequence[ ($position - $one_flanking_sequence_length) ... ($position + $one_flanking_sequence_length) ];
    }
    else {
        $self->warning_message("Length of wildtype sequence is shorter than desired peptide length of output. Skipping position $position");
    }

    return @wildtype_arr;
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

sub get_wildtype_sequence {
    my $self = shift;
    my @protein_arr = @_;

    if (scalar(@protein_arr) == 22) {
        return $protein_arr[21];
    }
    elsif (scalar(@protein_arr) == 25) {
        return $protein_arr[24];
    }
    else {
        confess $self->error_message("Protein array not as expected");
    }
}

sub print_fasta_entry {
    my $self = shift;
    my %params = @_;

    my $output_fh = Genome::Sys->open_file_for_appending($self->output_file);

    my ($fasta_defline, $fasta_sequence);
    my $identifier = $params{'designation'} . "." . $params{protein_arr}->[6] . "." . $params{protein_arr}->[15];
    $fasta_defline = ">$identifier";
    $fasta_sequence = join "", @{$params{arr}};

    print $output_fh "$fasta_defline\n";
    print $output_fh "$fasta_sequence\n";

    close($output_fh);
}

1;
