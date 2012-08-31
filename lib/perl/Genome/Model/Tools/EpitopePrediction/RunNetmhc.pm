package Genome::Model::Tools::EpitopePrediction::RunNetmhc;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::EpitopePrediction::RunNetmhc {
    is => ['Genome::Model::Tools::EpitopePrediction::Base'],
    has_input => [
        allele => {
        	is => 'Text',
            doc => 'Allele name to predict epitope prediction. For a list of available alleles, use: netMHC -A',
        },
        fasta_file => {
        	is => 'Text',
            doc => 'Input 21-mer Fasta file contains Wildtype(WT) and Mutant(MT) protein sequences',
        },
        output_file => {
            is => 'Text',
            doc => 'Output file containing raw output from Netmhc epitope prediction',
        },
        
        epitope_length => {
            is => 'Text',
            doc => 'Length of subpeptides to predict',
        },
    ],
};

sub help_brief {
    "FOR NETMHC : Runs NetMHC for MHC Class I epitope prediction

",
}



sub execute {
    my $self = shift;
    
    my ($temp_fh_name, $temp_name) = Genome::Sys->create_temp_file();
    
    my $netmhc_cmd = 'bsub -q techd -u jhundal@genome.wustl.edu -R \'select[mem >4000] rusage[mem=4000]\' -M 4000000 -N -oo'.
    				  $self->output_file.'_stdout'.
    				  ' /gsc/bin/netMHC -a '.
    				   $self->allele .
    				   ' -l '.$self->epitope_length.' '.
    				   $self->fasta_file.
    				   ' -x '.$self->output_file;

    Genome::Sys->shellcmd(
        cmd => $netmhc_cmd,
        input_files => [$self->fasta_file],
#        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    
    return 1;
}



1;
