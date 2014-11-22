package Genome::Model::Tools::UpdateTabCompletion;

use strict;
use warnings;

class Genome::Model::Tools::UpdateTabCompletion {
    is => 'Command::V2',
    doc => 'update the tab completion spec files (.opts)',
};

sub execute {
    my $self = shift;

    my @command_classes = ('Genome::Model::Tools', 'Genome::Command');
    for my $classname (@command_classes) {
        my $genome_completion = UR::Namespace::Command::CreateCompletionSpecFile->create(
            classname => $classname,
        );
        unless ($genome_completion->execute) {
            $self->error_message("Updating the $classname spec file did not complete succesfully!");
        }
    }

    $self->status_message("Remember to commit the updated .opts file(s)!");

    return 1;
}

1;
