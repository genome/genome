package Genome::Model::ReferenceSequence::Command::ListBuilds;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceSequence::Command::ListBuilds {
    is => 'Genome::Model::Command::BuildRelatedList',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Model::Build::ImportedReferenceSequence'
        },
        show => { default_value => 'id,name,model.subject.name,date_scheduled,status,run_by,data_directory' },
    ],
};

sub sub_command_sort_position { 1 }

sub help_synopsis {
    my $self = shift;
    my $syn = $self->SUPER::help_synopsis(@_);
    $syn .= <<EOS;

  # list all builds for the model "NCBI-human"
  genome model reference-sequence list-builds NCBI-human

  # list a specific build by name
  genome model reference-sequence list-builds name=NCBI-human-build36

  # or use standard filters
  genome model reference-sequence list-builds --filter name~UCSC%,status=Succeeded,subject_name="mouse" --show id,subject_name,data_directory

EOS
    return $syn;
}

1;

