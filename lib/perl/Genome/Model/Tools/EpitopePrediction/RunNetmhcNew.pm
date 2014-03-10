package Genome::Model::Tools::EpitopePrediction::RunNetmhcNew;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::EpitopePrediction::RunNetmhcNew {
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
        version => {
            is => 'Text',
            doc => 'NetMHC version to use',
            valid_values => ['3.0','3.4'],
            default_value => '3.4',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    "FOR NETMHC : Runs NetMHC for MHC Class I epitope prediction

",
}

sub execute {
    my $self = shift;

    my $netmhc_path;
    if ($self->version eq '3.0') {
        $netmhc_path = '/gsc/bin/netMHC';
    }
    elsif ($self->version eq '3.4') {
        $netmhc_path = '/gscmnt/sata141/techd/jhundal/netMHC/NetMHC3.4/ATTEMPT4/NetMHC/netMHC';
    }

    my $netmhc_cmd = join( ' ',
            'bsub -q techd -R \'select[mem >4000] rusage[mem=4000]\' -M 4000000 -N',
            '-oo ' . $self->output_file . '_stdout',
            $netmhc_path,
            '-a ' . $self->allele,
            '-l ' . $self->epitope_length,
            $self->fasta_file,
            '-x ' . $self->output_file,
    );

    Genome::Sys->shellcmd(
        cmd => $netmhc_cmd,
        input_files => [$self->fasta_file],
#        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );

    return 1;
}

1;
