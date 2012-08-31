package Genome::Sys::Command::Completion::Update;

use strict;
use warnings;

class Genome::Sys::Command::Completion::Update {
    is => 'Genome::Command::Base',
    doc => 'update the tab completion spec files (.opts)',
    has => [
        git_add => {
            is => 'Boolean',
            doc => 'git add the changed files after update',
            default => 0,
        },
        git_commit => {
            is => 'Boolean',
            doc => 'git commit the changed files after update',
            default => 0,
        },
    ],
};

sub help_detail {
    my $help_detail;

    $help_detail .= "Updates the tab completion spec files:\n";
    $help_detail .= " * Genome/Command.pm.opts\n";
    $help_detail .= " * Genome/Model/Tools.pm.opts";

    return $help_detail;
}

sub execute {
    my $self = shift;

    my @command_classes = ('Genome::Model::Tools', 'Genome::Command');
    for my $classname (@command_classes) {
        my $genome_completion = UR::Namespace::Command::Update::TabCompletionSpec->create(
            classname => $classname,
        );
        unless ($genome_completion->execute) {
            $self->error_message("Updating the $classname spec file did not complete succesfully!");
        }
    }

    my @files = qx(git status -s 2> /dev/null);
    chomp @files;

    my ($gmt_opts) = grep { $_ =~ /Tools\.pm\.opts/ } @files;
    $gmt_opts =~ s/^\s*M\s*// if $gmt_opts;

    my ($genome_opts) = grep { $_ =~ /Command\.pm\.opts/ } @files;
    $genome_opts =~ s/^\s*M\s*// if $genome_opts;

    my @opts_files;
    push @opts_files, $gmt_opts if $gmt_opts;
    push @opts_files, $genome_opts if $genome_opts;

    if (@opts_files) {
        if ($self->git_add) {
            my $rv = system('git add ' . join (' ', @opts_files));
            if ($rv == 0) {
                $self->status_message(join("\n\t", 'Added .opts file(s) to git:', @opts_files));
                $self->status_message("Remember to commit the .opts file(s).");
            }
        }
        elsif ($self->git_commit) {
            my $rv = system("git commit -m 'updated .opts file(s)' " . join(' ', @opts_files));
            if ($rv == 0) {
                $self->status_message(join("\n\t", 'Committed .opts file(s) to git:', @opts_files));
                $self->status_message("Remember to push the .opts file(s).");
            }
        }
        else {
            $self->status_message("Remember to commit and push the .opts file(s).");
        }
    }


    return 1;
}

1;
