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
            my ($wildtype_aa, $position, $mutant_aa) = ($1, $2 - 1, $3);
            my $wildtype_sequence;
            if (scalar(@protein_arr) == 22) {
                $wildtype_sequence = $protein_arr[21];
            }
            elsif (scalar(@protein_arr) == 25) {
                $wildtype_sequence = $protein_arr[24];
            }
            else {
                confess $self->error_message("Protein array not as expected");
            }
            my @arr_wildtype_sequence = split('',$wildtype_sequence);

            if ($wildtype_aa ne $arr_wildtype_sequence[$position]) {
                next;
            }
            else {
                my (@mutant_arr, @wildtype_arr);
                my $midpoint = ($self->peptide_sequence_length - 1) / 2;
                if ($position < $midpoint) {
                    @wildtype_arr = @mutant_arr = @arr_wildtype_sequence[ 0 ... ($self->peptide_sequence_length - 1) ];
                }
                elsif ($position > ($#arr_wildtype_sequence - $midpoint)) {
                    @wildtype_arr = @mutant_arr = @arr_wildtype_sequence[ ($#arr_wildtype_sequence - $self->peptide_sequence_length) ... $#arr_wildtype_sequence];
                }
                elsif (($position >= $midpoint) && ($position  <= ($#arr_wildtype_sequence - $midpoint))) {
                    @wildtype_arr = @mutant_arr = @arr_wildtype_sequence[ ($position - $midpoint) ... ($position + $midpoint) ];
                }
                else {
                    my $output_fh = Genome::Sys->open_file_for_appending($self->output_file);
                    print $output_fh "NULL"."\t".$position."\n";
                    close($output_fh);
                    next;
                }
                $mutant_arr[$midpoint]=$mutant_aa;

                for my $designation ('WT', 'MT') {
                    my $arr = $designation eq 'WT' ? \@wildtype_arr : \@mutant_arr;
                    $self->print_fasta_entry(
                        designation => $designation,
                        arr         => $arr,
                        protein_arr => \@protein_arr,
                    );
                }
            }
        }
    }

    close($input_fh);

    return 1;
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
