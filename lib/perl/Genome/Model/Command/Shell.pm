package Genome::Model::Command::Shell;
use strict;
use warnings;

use Genome;
class Genome::Model::Command::Shell {
    is => 'Command',
    doc => 'start a shell interface for genome model commands',
};

use Term::ReadLine;

sub sub_command_sort_position { 99 }

sub execute {
    my $self = shift; 

    my $term = Term::ReadLine->new('cool');
    #my $term = Term::ReadLine::Gnu->new('cool');

    #$term->bind_key('tab', \&tab_key);

    $SIG{'ALRM'} = sub { print "\nStill working...\n" };

    my $current_class = 'Genome::Command';
    my @class_stack = ();

    while(1) {
        alarm(0);

        my $prompt = $self->prompt_for_class($current_class);
        my $line = $term->readline($prompt);
        last unless defined $line;

        my @words = split(/\s+/, $line);

        if (@words == 1 and $words[0] eq '/') {
            $current_class = 'Genome::Command';
            @class_stack = ();
            next;
        } elsif (@words == 1 and $words[0] eq '..') {
            $current_class = pop @class_stack;
            $current_class = 'Genome::Command' unless $current_class;
            next;
        }

        alarm(2);

        my ($delegate_class, $params) = $current_class->resolve_class_and_params_for_argv(@words);

        unless ($delegate_class) {
            $current_class->usage_message($current_class->help_usage_complete_text);
            return 1;
        }

        if (!$params and ($current_class ne $delegate_class)) {
            push @class_stack, $current_class;
            $current_class = $delegate_class;
            next;
        }
        if (!$params or $params->{help}) {
            $delegate_class->usage_message($delegate_class->help_usage_complete_text,"\n");
            next;
        }

        if ($delegate_class eq __PACKAGE__ ) {
            $self->error_message("We're not letting you start a shell inside a shell");
            next;
        }

        my $command = $delegate_class->create(%$params);
        unless ($command) {
            print STDERR "Error\n";
            next;
        }

        my $rv = $command->execute($params);
        if ($command->__errors__) {
            $command->delete;
            print STDERR "Invalid!?\n";
            next;
        }

        if ($rv) {
            print "SUCCESS\n";
            UR::Context->commit();
        } else {
            print "FAILED\n";
            UR::Context->rollback();
        }
    }
    return 1;
}

sub prompt_for_class {
    my($self,$path) = @_;

    my @parts = map { lcfirst }
                split(/::/,$path);
    return join(' ',@parts) . '> ';
}

1;
