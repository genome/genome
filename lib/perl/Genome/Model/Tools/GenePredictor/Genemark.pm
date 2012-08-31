package Genome::Model::Tools::GenePredictor::Genemark;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use IO::Dir;

# TODO Remove this dependency
use BAP::Job::Genemark;

class Genome::Model::Tools::GenePredictor::Genemark {
    is => ['Genome::Model::Tools::GenePredictor'],
    has => [
        gc_percent => { 
            is => 'Number', 
            doc => 'GC content',
            is_input => 1, 
        },
    ],
    has_optional => [
        model_file => { 
            is => 'FilePath',
            doc => 'Genemark model file',
        },
    ],
};

sub help_brief {
    "Write a set of fasta files for an assembly";
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

    my $seqio = Bio::SeqIO->new(-file => $self->fasta_file(), -format => 'Fasta');

    my $seq = $seqio->next_seq();

    my $gc_percent = sprintf("%.0f", $self->gc_percent());
    my $gc_model = $self->_select_model($gc_percent);

    $self->model_file($gc_model);
    
    ##FIXME: The last arg is the job_id, which is hardcoded here in a rather lame fashion.
    # TODO Remove this rather lame job, just execute the command here without using this extra
    # layer of indirection
    my $legacy_job = BAP::Job::Genemark->new(
                                             $seq,
                                             $gc_model,
                                             2112,
                                        );

    $legacy_job->execute();

    my @features = $legacy_job->seq()->get_SeqFeatures();
    
    $self->{bio_seq_feature} = \@features;
           
    return 1;
    
}

sub _select_model {

    my $self       = shift;
    my $gc_percent = shift;
    

    ##FIXME: This should probably not be hardcoded, at least not here
    my $model_dir = $ENV{GENOME_SW} . '/genemark.hmm/installed/modeldir';

    ##FIXME: These thresholds should not be hardcoded, either
    if ($gc_percent < 30) {
        $gc_percent = 30;
    }
    elsif ($gc_percent > 70) {
        $gc_percent = 70;
    }
    
    unless (-e $model_dir) {
        die "model directory does not seem to exist: $model_dir";
    }

    my $dh = IO::Dir->new($model_dir);

    my @model_files = $dh->read();

    my ($model_file) = grep { $_ =~ /heu_11_$gc_percent\.mod/ } @model_files;

    unless (defined($model_file)) {
        die "could not locate model file for gc content of $gc_percent percent";
    }

    return "$model_dir/$model_file";
    
}

1;
