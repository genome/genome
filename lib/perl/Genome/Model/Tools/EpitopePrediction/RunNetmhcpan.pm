package Genome::Model::Tools::EpitopePrediction::RunNetmhcpan;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::EpitopePrediction::RunNetmhcpan {
    is => ['Genome::Model::Tools::EpitopePrediction::Base'],
    has_input => [
        allele => {
        	is => 'Text',
            doc => 'Allele name to predict epitope prediction. For a list of available alleles, use: netMHCpan -listMHC',
        },
        fasta_file => {
        	is => 'Text',
            doc => 'Input 21-mer Fasta file containing Wildtype(WT) and Mutant(MT) protein sequences',
        },
        output_file => {
            is => 'Text',
            doc => 'Output file containing raw output from NetmhcPan epitope prediction',
        },
        
        epitope_length => {
            is => 'Text',
            doc => 'Length of subpeptides to predict',
        },
    ],
};

sub help_brief {
    "FOR NETMHCPAN : Runs NetMHCPAN for MHC Class I epitope prediction

",
}



sub execute {
    my $self = shift;
    
    #netMHCpan -a HLA-A03:01 -f WHIM32_21_ns.fa -xlsfile WHIM32_21_pan.xls -xls > pan_netmhc
    
    my $netmhcpan_cmd = 'bsub -q techd -R \'select[mem >8000] rusage[mem=8000]\' -M 8000000 -N '.
    				  
    				  '\'/gsc/bin/netMHCpan -a '.
    				   $self->allele .
    				   ' -l '.$self->epitope_length.' '.
    				   ' -f '.$self->fasta_file." > ".$self->output_file.'\'' ;
	
print $netmhcpan_cmd ;
   			   

    Genome::Sys->shellcmd(
       						 cmd => $netmhcpan_cmd,
  						  #    input_files => [$self->fasta_file],
    #    output_files => [$self->output_file],
    #   skip_if_output_is_present => 0,
    );
    
    return 1;
}



1;
