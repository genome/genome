package Genome::Model::Tools::Annotate::GenerateTranscriptCodingSequenceFile;

use strict;
use warnings;

use Genome;

#TODO: change the UR object cache limits so they are higher

class Genome::Model::Tools::Annotate::GenerateTranscriptCodingSequenceFile{
    is => 'Genome::Model::Tools::Annotate',
    has => [
        reference_transcripts => {
            is => 'String',
            is_input => 1, 
            is_optional => 0,
            doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/0")',
        },
    ],
};

sub help_synopsis{
    return <<EOS
gmt annotate generate-transcript-coding-sequence-file --reference_transcripts NCBI-human.combined-annotation/54_36p
EOS
}

#TODO: Do something helpful
sub help_detail{
    return "Generates a transcript coding sequence file from a model/build string";
}

sub execute{
    my $self = shift;

    my ($model_name, $build_version) = split("/", $self->reference_transcripts);
    my $model = Genome::Model->get(name => $model_name);
    die "Could not get model $model_name" unless $model;
    my $build = $model->build_by_version($build_version);
    die "Could not get imported annotation build version $build_version" unless $build;

    my @builds = $build->from_builds;
    if (@builds) {
        for my $b (@builds) {
            $self->create_sequence_file($b);
        }
    }
    else {
        $build->create_sequence_file($build);
    }

    return 1;
}

sub create_sequence_file {
    my ($self, $build) = @_;

    my $iterator = $build->transcript_iterator;
    die "Could not get iterator" unless $iterator;
    while ( my $transcript = $iterator->next) {
        my $seq = $transcript->cds_full_nucleotide_sequence;
        next unless defined $seq;
        my $id = $transcript->id;
        my $obj = Genome::TranscriptCodingSequence->create(
            sequence => $seq,
            transcript_id => $id,
            data_directory => $transcript->data_directory,
        );
    }
}
1;
