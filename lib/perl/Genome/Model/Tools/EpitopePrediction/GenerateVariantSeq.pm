package Genome::Model::Tools::EpitopePrediction::GenerateVariantSeq;

use strict;
use warnings;

use Genome;
use Workflow;

my $DEFAULT_TRV_TYPE = 'missense';


class Genome::Model::Tools::EpitopePrediction::GenerateVariantSeq {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input file is a tab-separated (TSV) output file from gmt annotate variant-protein. For more info, gmt annotate variant-protein --help',
        },
        output_file => {
            is => 'Text',
            doc => 'The output FASTA file to write 21mer sequences for wildtype(WT) and mutant(MT) proteins',
        },
        trv_type => {
            is => 'Text',
            is_optional => 1,
            doc => 'The type of mutation you want to output eg missense,nonsense',
            # the current code only works on missense. Will need furthur development for other trv_types.
            default_value => $DEFAULT_TRV_TYPE,
        },
        length => {
            is => 'Number',
            doc => 'The length of the peptide sequences',
            valid_values => [17, 21, 31],
        },
    ],
};


sub help_brief {
    "FOR NETMHC : Outputs a FASTA file for  for wildtype(WT) and mutant(MT) proteins 21-mer sequences for MHC Class I epitope prediction",
}



sub execute {
    my $self = shift;
    #my $tmp_dir = Genome::Sys->create_temp_directory();

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
                #TO DO :print $output_fh $protein_arr[0]."\t".$protein_arr[1]."\t".$protein_arr[2]."\t".$protein_arr[6]."\t".$1."\t".$2."\t".$3."\t".$protein_arr[11]."\t".$arr_wildtype_sequence[$position]."\n";
            }
            else {
                my (@mutant_arr, @wildtype_arr);
                my $midpoint = ($self->length - 1) / 2;
                if ($position < $midpoint) {
                    @wildtype_arr = @mutant_arr = @arr_wildtype_sequence[ 0 ... ($self->length - 1) ];
                }
                elsif ($position > ($#arr_wildtype_sequence - $midpoint)) {
                    @wildtype_arr = @mutant_arr = @arr_wildtype_sequence[ ($#arr_wildtype_sequence - $self->length) ... $#arr_wildtype_sequence];
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
