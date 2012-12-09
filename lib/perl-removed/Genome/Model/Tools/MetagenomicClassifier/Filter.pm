package Genome::Model::Tools::MetagenomicClassifier::Filter;
use strict;
use warnings;

require Bio::SeqIO;
use Genome;
require Genome::Utility::MetagenomicClassifier::Rdp;
require Genome::Utility::MetagenomicClassifier::ChimeraClassifier;
require Genome::Utility::MetagenomicClassifier::Weka;
require IO::File;


class Genome::Model::Tools::MetagenomicClassifier::Filter {
    is => 'Command',
    has => [ 
        input_file => {
            type => 'String',
            doc => "path to fasta file"
        },
    ],
    has_optional => [
        training_path => {
            type => 'String',
            doc => 'path to training set directory',
        },
        training_set => {
            type => 'String',
            doc => 'name of training set in default location (broad)',
        },
        
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
    
    my $clean_detail_writer = Genome::Utility::MetagenomicClassifier::SequenceClassification::Writer->create(
        output => $self->input_file.".clean.detail",
        format => 'hmp_all_ranks',
    );

    my $chimera_detail_writer =  Genome::Utility::MetagenomicClassifier::ChimeraClassification::Writer->create(
        verbose => 1,
        output => $self->input_file.".chimera.detail",
    );

    my $unclassified_detail_writer =  Genome::Utility::MetagenomicClassifier::ChimeraClassification::Writer->create(
        verbose => 1,
        output => $self->input_file.".unclassified.detail",
    );

    
    my $weka = Genome::Utility::MetagenomicClassifier::Weka->create();
    #< IN >#
    my $bioseq_in = Bio::SeqIO->new(
        -format => 'fasta',
        -file => $self->input_file,
    ) or return;

    #< OUT >#
    my $clean_writer = Bio::SeqIO->new(
        -format => 'fasta',
        -file => ">".$self->input_file.".clean.fasta",
    );

    my $chimera_writer = Bio::SeqIO->new(
        -format => 'fasta',
        -file => ">".$self->input_file.".chimera.fasta",
    );

    my $unclassified_writer = Bio::SeqIO->new(
        -format => 'fasta',
        -file => ">".$self->input_file.".unclassified.fasta",
    );

    while ( my $seq = $bioseq_in->next_seq ) {
        my $classification = $classifier->classify($seq);
        if ($classification) {
            my $class = $weka->classify($classification->get_profile);
            if ($class eq 'clean') {
                $clean_writer->write_seq($seq);
                $clean_detail_writer->write_one($classification->classification);
            }
            elsif ($class eq 'chimera') {
                $chimera_writer->write_seq($seq);
                $chimera_detail_writer->write_one($classification);
            }
            else {
                $unclassified_writer->write_seq($seq);
                $unclassified_detail_writer->write_one($classification);
            } 
        }

    }

    return 1;
}

#< HELP >#
sub help_brief {
    "interesting sequence filter",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools metagenomic-classifier filter    
EOS
}

1;

#$HeadURL$
#$Id$
