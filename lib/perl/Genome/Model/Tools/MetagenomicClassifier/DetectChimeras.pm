package Genome::Model::Tools::MetagenomicClassifier::DetectChimeras;
use strict;
use warnings;

use Bio::SeqIO;
use Genome::Utility::MetagenomicClassifier::ChimeraClassifier;
use IO::File;

class Genome::Model::Tools::MetagenomicClassifier::DetectChimeras {
    is => 'Command',
    has => [ 
        input_file => {
            type => 'String',
            doc => "path to fasta file"
        },
    ],
    has_optional => [
        output_file => { 
            type => 'String',
            is_optional => 1, ###
            doc => "path to output file.  Defaults to STDOUT"
        },
        training_path => {
            type => 'String',
            doc => 'path to training set directory',
        },
        training_set => {
            type => 'String',
            doc => 'name of training set in default location (broad)',
        },
        verbose => {
            type => 'BOOL',
            default => 0,
            doc => 'output the full classification for each probe',
        },
        arff => {
            type => 'BOOL',
            default => 0,
            doc => 'output arff header',
        }
    ],
};

sub execute {
    my $self = shift;
    
    my $training_path = '/gsc/scripts/share/rdp/';
    if ($self->training_path) {
        $training_path = $self->training_path.'/';
    }
    elsif ($self->training_set) {
        $training_path .= $self->training_set.'/';
    }

    my $rdp_classifier = Genome::Utility::MetagenomicClassifier::Rdp->new(
        training_path => $training_path,
    );
    #< CLASSIFER >#
    my $classifier = Genome::Utility::MetagenomicClassifier::ChimeraClassifier->create(
        classifier => $rdp_classifier,
    ) or return;
    
    #< IN >#
    my $bioseq_in = Bio::SeqIO->new(
        -format => 'fasta',
        -file => $self->input_file,
    ) or return;

    #< OUT >#
    my $writer = Genome::Utility::MetagenomicClassifier::ChimeraClassification::Writer->create(
        output => $self->output_file,
        verbose => $self->verbose,
        arff => $self->arff,
    )
        or return;

    if ($self->arff) {
        $writer->write_arff_header;
    }
    while ( my $seq = $bioseq_in->next_seq ) {
        my $classification = $classifier->classify($seq);
        $writer->write_one($classification) if $classification;
    }

    return 1;
}

#< HELP >#
sub help_brief {
    "chimera detector",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools metagenomic-classifier detect-chimeras    
EOS
}

1;

#$HeadURL$
#$Id$
