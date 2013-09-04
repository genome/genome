package Genome::Model::Command::Diff;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::Command::Diff {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model',
            require_user_verify => 0,
            shell_args_position => 1,
            is_many => 1,
            doc => 'Models that should have their builds compared',
        },
        first_revision => {
            is => 'Text',
            doc => 'Path to revision that one build was run on (eg, /gsc/scripts/opt/genome/current/stable)',
        },
        second_revision => {
            is => 'Text',
            doc => 'Path to revision that other build was run on'
        },
    ],
};

sub hudson_build_from_revision {
    my ($self, $revision) = @_;
    if ($revision =~ /(genome-\d+)/) {
        return $1;
    }
    return $revision;
}

# Check that the revision ends with lib/perl and add it if its not there
sub check_and_fix_revision {
    my ($self, $revision) = @_;
    $revision .= '/' unless $revision =~ /\/$/;
    $revision .= 'lib/perl' unless $revision =~ 'lib/perl';
    return $revision;
}

sub execute { 
    my $self = shift;

    my $first_revision = $self->first_revision;
    $first_revision = readlink $first_revision if -l $first_revision;
    $first_revision = $self->check_and_fix_revision($first_revision);
    confess "Revision not found at $first_revision!" unless -d $first_revision;

    my $second_revision = $self->second_revision;
    $second_revision = readlink $second_revision if -l $second_revision;
    $second_revision = $self->check_and_fix_revision($second_revision);
    confess "No revision found at $second_revision!" unless -d $second_revision;

    if ($first_revision eq $second_revision) {
        confess "Comparing builds from $first_revision and $second_revision... these are the same, no point in comparing.";
    }

    $self->status_message("Comparing builds from revisions $first_revision and $second_revision");

    # If comparing hudson nightly builds, the software revision of the build with be /gsc/scripts/opt/passed-unit-tests*,
    # which no longer exists at this point because a successful model test results in the snapshot directory being moved
    # from passed-unit-tests to passed-model-tests. Using only the genome-### portion is good enough.
    my $fixed_first_revision = $self->hudson_build_from_revision($first_revision);
    my $fixed_second_revision = $self->hudson_build_from_revision($second_revision);
    for my $model ($self->models) {
        my $model_id = $model->genome_model_id;
        my $type = $model->class;
        $type =~ s/Genome::Model:://;
        $type =~ s/:://g;
        if ($type =~ /Convergence/) {
            # Talk to Tom for details... basically, there's no expectation that the output be
            # the same between builds, so diffing the output at all is pointless.
            $self->status_message("This is a convergence model. Skipping...");
            next;
        }
                                        
        my $type_string = Genome::Utility::Text::camel_case_to_string($type, '_');
        $self->status_message("\nWorking on model " . $model->name . " ($model_id), type $type_string");

        # Find builds for each given revision
        my ($first_build, $second_build) = $self->find_builds_for_revisions($model, $fixed_first_revision, $fixed_second_revision);
        unless ($first_build and $second_build) {
            my $msg = "BUILD NOT FOUND $type_string $model_id: Could not find build for model $model_id using revision:";
            $msg .= ' ' . $fixed_first_revision unless $first_build;
            $msg .= ' ' . $fixed_second_revision unless $second_build;
            $self->warning_message($msg);
            next;
        }

        $self->status_message("Comparing build " . $first_build->build_id . " using revision $first_revision with build " . 
            $second_build->build_id . " using revision $second_revision from model $model_id");
        $self->status_message("Build " . $first_build->build_id . " has data directory " . $first_build->data_directory .
            " and build " . $second_build->build_id . " has data directory " . $second_build->data_directory);

        my %diffs = $first_build->compare_output($second_build->build_id);
        unless (%diffs) {
            $self->status_message("All files and metrics diffed cleanly!");
        }
        else {
            my $diff_string = "DIFFERENCES FOUND $type_string $model_id\n";
            for my $file (sort keys %diffs) {
                my $reason = $diffs{$file};
                $diff_string .= "  File: $file\n  Reason: $reason\n";
            }
            $self->status_message($diff_string);
        }
    }

    return 1;
}

sub find_builds_for_revisions {
    my ($self, $model, $first_revision, $second_revision) = @_;

    my @builds = sort { $b->date_scheduled <=> $a->date_scheduled } grep { $_->status eq 'Succeeded' } $model->builds;
    my $user   = Genome::Sys->username;
    my ($first_build, $second_build);

    for my $build (@builds) {
        if ($user eq 'apipe-tester') {
            next unless $build->run_by eq $user;
        }
        last if $first_build and $second_build;
        my $build_revision = $build->software_revision;
        next unless defined $build_revision;
        if (not $first_build and $build_revision =~ /$first_revision/) {
            $first_build = $build;
        }
        elsif (not $second_build and $build_revision =~ /$second_revision/) {
            $second_build = $build;
        }
    }

    return ($first_build, $second_build);
}

1;  

