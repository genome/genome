#$Id$

package Genome::Model::GenePrediction::Command::Pap::Hmmpfam;

use strict;
use warnings;

use English;
use File::Basename;
use File::chdir;
use File::Temp;
use File::Slurp;
use IO::File;
use IPC::Run;
use Carp;

class Genome::Model::GenePrediction::Command::Pap::Hmmpfam {
    is  => ['Command::V1'],
    has => [
        hmmdirpath => {
            is  => 'SCALAR',
            doc => 'path to hmm output files',
            is_input => 1,
        },
        hmm_database => {
            is  => 'SCALAR',
            doc => 'hmmpfam database file',
            is_input => 1,
        },
        fasta_dir => {
            is  => 'SCALAR',
            doc => 'directory containing fasta files',
            is_input => 1,
        },
        sequence_names => {
            is  => 'ARRAY',
            doc => 'a list of sequence names for running hmmpfam on',
            is_input => 1,
        },
        success => {
            is          => 'SCALAR',
            doc         => 'success flag',
            is_optional => 1,
            is_output => 1,
        },
        lsf_queue => { 
            is_param => 1, 
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
        lsf_resource => { 
            is_param => 1,
            default_value => '-R "rusage[tmp=100]" ',
        },
    ],
};


sub sub_command_sort_position {10}

sub help_brief
{
    "Run run the BER hmmpfam step";
}

sub help_synopsis
{
    return <<"EOS"
EOS
}

sub help_detail
{
    return <<"EOS"
Need documenation here.
EOS
}

sub execute
{

    my $self = shift;

    {

#        local $CWD = $self->workdir();
        foreach my $seqname ( @{ $self->sequence_names } )
        {
            my $file = $self->fasta_dir . "/" . $seqname;

            my @hmmpfam_command
                = ( 'hmmpfam', '--cpu', '1', $self->hmm_database, $file, );

            my ( $hmmpfam_stdout, $hmmpfam_stderr );
            my $hmmpfam_outfile
                = $self->hmmdirpath . "/" . $seqname . ".hmmpfam";
            IPC::Run::run(
                \@hmmpfam_command, 
                '<',  
                \undef, 
                '>', \$hmmpfam_stdout,  
                '2>', \$hmmpfam_stderr,
                )
                || croak "ber hmmpfam failed: $hmmpfam_stderr : $CHILD_ERROR";
            write_file( $hmmpfam_outfile, $hmmpfam_stdout );
        }
    }
    $self->success(1);
    return 1;

}

1;
