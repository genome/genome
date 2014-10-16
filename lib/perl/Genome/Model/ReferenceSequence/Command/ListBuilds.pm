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
  # given model 123 named "mymodel" has builds 456 and 789:

  # list the first build
  genome model build list 456

  # list the second build
  genome model build list 789

  # list all builds for the model by model ID
  genome model build list 123

  # list all builds for the model by model name
  genome model build list mymodel

  # or use standard filters
  genome model build list --filter status=Abandoned,data_directory~/gscmnt/839% --show id,subject_name,data_directory
EOS
    return $syn;
}

1;

