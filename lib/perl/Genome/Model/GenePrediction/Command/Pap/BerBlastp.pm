#$Id$

package Genome::Model::GenePrediction::Command::Pap::BerBlastp;

use strict;
use warnings;

#use Workflow;

use Bio::Annotation::DBLink;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;

use English;
use File::Basename;
use File::Spec;
use File::Temp;
use File::Slurp;
use IPC::Run;
use Carp;


class Genome::Model::GenePrediction::Command::Pap::BerBlastp {
    is  => ['Command::V1'],
    has => [

        sequence_names => { is => 'ARRAY',
                            doc => 'array of sequence names',
                            is_input => 1,
                          },
        fastadir => { is => 'SCALAR',
                      doc => 'location of fasta files',
                      is_input => 1,
                    },
        blastp_query => { is => 'SCALAR',
                          doc => 'queryfile for blastp',
                          is_input => 1,
                        },
        berdirpath => { is => 'SCALAR',
                        doc => 'directory where blastp output lands',
                        is_input => 1,
                      },
        success => { is => 'SCALAR',
                     doc => 'successful execution flag',
                     is_optional => 1,
                     is_output => 1,
                   },
        lsf_queue => { is_param => 1,
                       default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
        lsf_resource => { is_param => 1,
                          default_value => 'rusage[tmp=100]',
        },
    ],
};

#operation Genome::Model::GenePrediction::Command::Pap::BerBlastp {
#    input        => [ 'sequence_names', 'fastadir','blastp_query', 'berdirpath' ],
#    output       => [ 'success'],
#    lsf_queue    => 'long',
#    lsf_resource => 'rusage[tmp=100]'
#};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Run BER blastp step";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {

    my $self = shift;

    my $blp_err_path = $self->berdirpath;
    $blp_err_path =~ s/ber/bsubERRfiles/;
    $self->status_message('getting ready to run this stuff');
    unless(-d $blp_err_path)
    {
        warn "creating directory $blp_err_path for error files";
        mkdir($blp_err_path);
    }

    $self->status_message('starting run of blastp...');
    foreach my $seqname (@{$self->sequence_names})
    {

        my ($blastp_out, $blastp_err);
        my $fasta_file = $self->fastadir() . "/".$seqname;
    
        ## If 'blastp' invokes anything but WU-BLAST, stuff will probably
        ## go seriously foul in archive_result below
        ## it looks like there were not any extra options or params
        ## added to the blastp query for BER...
        my @blastp_command = (
                              'blastp',
                              $self->blastp_query,
                              $fasta_file,
                             );

        IPC::Run::run(
                      \@blastp_command,
                      \undef,
                      '>',
                      \$blastp_out,
                      '2>',
                      \$blastp_err,
                  ) || croak "ber blastp failed: $CHILD_ERROR";
    
        my $blast_output_file = $self->berdirpath ."/".$seqname.".nr";
        write_file($blast_output_file,$blastp_out);
        write_file($blp_err_path."/bsub.err.blp.".$seqname,
                   $blastp_err);
    }    
    $self->status_message('blastp done.');
    $self->success(1);
    return 1;

}


1;
