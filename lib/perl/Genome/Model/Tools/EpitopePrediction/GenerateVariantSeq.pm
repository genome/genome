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
    ],
};


sub help_brief {
    "FOR NETMHC : Outputs a FASTA file for  for wildtype(WT) and mutant(MT) proteins 21-mer sequences for MHC Class I epitope prediction",
}



sub execute {
    my $self = shift;
    #my $tmp_dir = Genome::Sys->create_temp_directory();

    my $input_fh = Genome::Sys->open_file_for_reading($self->input_file);
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);

    while (my $line = $input_fh->getline) {
        chomp $line;
        $line =~ s/[*]$//g;
        my @protein_arr =  split(/\t/, $line);
        if ( $protein_arr[15] =~ /^p.([A-Z])(\d+)([A-Z])/ && $protein_arr[13] eq $self->trv_type) {
            my $wildtype_aa = $1;
            my $position = ($2 - 1);
            my $mutant_aa = $3;
            my $wildtype_sequence = $protein_arr[21];
            my @arr_wildtype_sequence = split('',$wildtype_sequence);

            if ($wildtype_aa ne $arr_wildtype_sequence[$position]) {
                next;
                #TO DO :print $output_fh $protein_arr[0]."\t".$protein_arr[1]."\t".$protein_arr[2]."\t".$protein_arr[6]."\t".$1."\t".$2."\t".$3."\t".$protein_arr[11]."\t".$arr_wildtype_sequence[$position]."\n";
            }
            else {
                my @mutant_arr;
                my @wildtype_arr;
                if ($position < 8) {
                    @wildtype_arr = @arr_wildtype_sequence[ 0 ... 16];
                    $arr_wildtype_sequence[$position]=$mutant_aa;
                    @mutant_arr = @arr_wildtype_sequence[ 0 ... 16];
                    print $output_fh ">WT.".$protein_arr[6].".".$protein_arr[15]."\n";
                    print $output_fh ( join "", @wildtype_arr);
                    print $output_fh "\n";
                    print $output_fh ">MT.".$protein_arr[6].".".$protein_arr[15]."\n";
                    print $output_fh ( join "", @mutant_arr);
                    print $output_fh "\n";
                }
                elsif ($position > ($#arr_wildtype_sequence -8)) {
                    @wildtype_arr = @arr_wildtype_sequence[ $#arr_wildtype_sequence -17 ... $#arr_wildtype_sequence];
                    $arr_wildtype_sequence[$position]=$mutant_aa;
                    @mutant_arr = @arr_wildtype_sequence[ $#arr_wildtype_sequence -17 ... $#arr_wildtype_sequence];
                    print $output_fh ">WT.".$protein_arr[6].".".$protein_arr[15]."\n";
                    print $output_fh ( join "", @wildtype_arr);
                    print $output_fh "\n";
                    print $output_fh ">MT.".$protein_arr[6].".".$protein_arr[15]."\n";
                    print $output_fh ( join "", @mutant_arr);
                    print $output_fh "\n";
                }
                elsif (($position >= 8) && ($position  <= ($#arr_wildtype_sequence -8))) {
                    @wildtype_arr = @arr_wildtype_sequence[ $position-8 ... $position+8];
                    $arr_wildtype_sequence[$position]=$mutant_aa;
                    @mutant_arr = @arr_wildtype_sequence[ $position-8 ... $position+8];
                    print $output_fh ">WT.".$protein_arr[6].".".$protein_arr[15]."\n";
                    print $output_fh ( join "", @wildtype_arr);
                    print $output_fh "\n";
                    print $output_fh ">MT.".$protein_arr[6].".".$protein_arr[15]."\n";
                    print $output_fh ( join "", @mutant_arr);
                    print $output_fh "\n";
                }
                else {
                    print $output_fh "NULL"."\t".$position."\n";
                }
            }
        }
    }

    close($output_fh);
    close($input_fh);

    return 1;
}

1;
