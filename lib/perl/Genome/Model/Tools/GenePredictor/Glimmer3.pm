package Genome::Model::Tools::GenePredictor::Glimmer3;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;

# TODO Remove this dependency
use BAP::Job::Glimmer;

class Genome::Model::Tools::GenePredictor::Glimmer3 {
    is => ['Genome::Model::Tools::GenePredictor'],
    has => [
        model_file => { 
            is => 'FilePath', 
            doc => 'absolute path to the model file for this fasta',
            is_input => 1, 
        },
        pwm_file => { 
            is => 'FilePath' , 
            doc => 'absolute path to the pwm file for this fasta', 
            is_input => 1, 
        },
    ],
};

sub help_brief {
    "Write a set of fasta files for an assembly";
}

sub execute {
    my $self = shift;
    my $seqio = Bio::SeqIO->new(-file => $self->fasta_file(), -format => 'Fasta');
    # TODO It'd be nice (and easy) to support multi-sequence fasta files
    my $seq = $seqio->next_seq();

    ##FIXME: The last two args are the circular dna flag and the
    ##       job_id, which are hardcoded here in a rather lame fashion.
    # TODO Remove dependency on lame jobs
    my $legacy_job = BAP::Job::Glimmer->new(
                                            'glimmer3',
                                            $seq,
                                            $self->model_file(),
                                            $self->pwm_file(),
                                            0,
                                            2112,
                                        );

    $legacy_job->execute();
    my @features = $legacy_job->seq()->get_SeqFeatures();
    $self->{bio_seq_feature} = \@features;
           
    return 1;
    
}

1;
