package Genome::Model::Tools::EpitopePrediction::RunNetmhc;

use strict;
use warnings;

use Genome;
use IPC::Run qw(run);

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
        output_directory => {
            is => 'Text',
            doc => 'Location of the output',
        },
        epitope_length => {
            is => 'Text',
            doc => 'Length of subpeptides to predict',
        },
        netmhc_version => {
            is => 'Text',
            doc => 'NetMHC version to use',
            valid_values => ['3.0','3.4'],
            default_value => '3.4',
            is_optional => 1,
        },
        sample_name => {
            is => 'Text',
            doc => 'The sample name to use in the name of the output file',
        },
    ],
    has_output => [
        output_file => {
            is => 'Text',
            doc => 'Output file containing raw output from Netmhc epitope prediction',
            calculate_from => ['output_directory', 'sample_name', 'allele', 'epitope_length'],
            calculate => q| return File::Spec->join($output_directory, "$sample_name.$allele.$epitope_length.netmhc.xls"); |,
        },
    ],
    has => [
        stdout_file => {
            is => 'Text',
            doc => 'Stdout file from Netmhc epitope prediction',
            calculate_from => ['output_directory', 'allele', 'epitope_length'],
            calculate => q| return File::Spec->join($output_directory, "$allele.$epitope_length.netmhc.stdout"); |,
        },

    ],
};

sub help_brief {
    "FOR NETMHC : Runs NetMHC for MHC Class I epitope prediction

",
}

sub execute {
    my $self = shift;

    my $netmhc_path = $self->netmhc_path;

    my @netmhc_cmd = (
            $netmhc_path,
            '-a', $self->allele,
            '-l', $self->epitope_length,
            $self->fasta_file,
            '-x', $self->output_file,
    );

    run(\@netmhc_cmd, ">", $self->stdout_file);

    $self->validate_output();

    return 1;
}

sub netmhc_path {
    my $self = shift;

    my %netmhc_path_of_version = (
        '3.0' => '/gsc/bin/netMHC',
        '3.4' => '/gscmnt/sata141/techd/jhundal/netMHC/NetMHC3.4/ATTEMPT4/NetMHC/netMHC',
    );

    return $netmhc_path_of_version{$self->netmhc_version};
}

sub validate_output {
    my $self = shift;

    unless (-s $self->output_file) {
        die $self->error_message("Output file %s does not exist or has no size after processing", $self->output_file);
    };

    unless (-s $self->stdout_file) {
        die $self->error_message("Stdout file %s does not exist or has no size after processing", $self->stdout_file);
    };
}

1;
