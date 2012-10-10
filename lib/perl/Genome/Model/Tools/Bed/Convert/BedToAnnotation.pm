package Genome::Model::Tools::Bed::Convert::BedToAnnotation;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Bed::Convert::BedToAnnotation{
    is => ['Genome::Model::Tools::Bed'],
    has => [
        snv_file => {
            is => 'Path',
            is_input => 1,
            is_optional => 1,
            doc => 'SNV file in bed format',
        },
        indel_file => {
            is => 'Path',
            is_input => 1,
            is_optional => 1,
            doc => 'indel file in bed format',
        },
        output => {
            is => 'Path',
            is_input => 1,
            is_output => 1,
            is_optional => 0,
            doc => 'File where the converted variants will be written',
        },
        annotator_version => {
          doc => 'Which version of the annotator\'s BedToAnnotation will be called',
          is => 'Text',
          valid_values => Genome::Model::Tools::Annotate::TranscriptVariants->__meta__->property(property_name => 'use_version')->valid_values,
          default_value => '1',
          is_input => 1,
          is_output => 0,
        },
    ],
};

sub execute{
    my $self = shift;

    my $annotation_version = "Version" . $self->annotator_version;
    my $command_class = "Genome::Model::Tools::Annotate::TranscriptVariants::${annotation_version}::BedToAnnotation";

    my $cmd = $command_class->create(
      snv_file   => $self->snv_file,
      indel_file => $self->indel_file,
      output     => $self->output,
    );

    $cmd->execute() || die("Failed to execute command...");

    return 1;
}

1;
