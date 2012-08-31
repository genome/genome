package Genome::Model::Tools::CompareBuilds;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::Tools::CompareBuilds {
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
            is => 'DirectoryPath',
            doc => 'Software revision that one of the builds should have been run on',
        },
    ],
    has_optional => [
        second_revision => {
            is => 'DirectoryPath',
            doc => 'Software revision that the other build should be run on, defaults to current/pipeline',
        },
    ],
};

sub hudson_build_from_software_revision {
    my ($self, $path) = @_;
    confess "Need a path to determine hudson build ID!" unless defined $path;
    $path =~ /(genome-\d+)/;
    return $1;
}

sub execute { 
    my $self = shift;

    # Can't use full path for snapshot here due to how hudson runs nightly builds... the software revision
    # for those builds is /gsc/scripts/opt/passed-unit-tests/genome-#/lib/perl, but the snapshot is then
    # moved to /gsc/scripts/opt/passed-model-tests/genome-#/lib/perl if builds are successful. Since the
    # genome-# part is unique and unchanging, this is sufficient to identify the correct revision.
    my $first_revision = $self->first_revision;
    confess "No revision found at $first_revision" unless -d $first_revision;
    $first_revision = $self->hudson_build_from_software_revision($first_revision);
   
    my $second_revision = $self->second_revision;
    unless (defined $second_revision) {
        my $stable = readlink "/gsc/scripts/opt/genome/current/pipeline";
        $second_revision = "/gsc/scripts/opt/$stable";
    }
    confess "No revision found at $second_revision" unless -d $second_revision;
    $second_revision = $self->hudson_build_from_software_revision($second_revision);

    if ($first_revision eq $second_revision) {
        confess "Comparing builds from $first_revision and $second_revision... these are the same, no point in comparing.";
    }

    $self->status_message("Comparing builds from revisions $first_revision and $second_revision");

    for my $model ($self->models) {
        my $model_id = $model->genome_model_id;
        my $type = $model->class;
        $type =~ s/Genome::Model:://;
        next if $type =~ /Convergence/; # Talk to Tom for details... basically, there's no expectation that the output be
                                        # the same between builds, so diffing the output at all is pointless.
    
        $self->status_message("\nWorking on model $model_id, type $type");

        my ($first_build, $second_build);
        my @builds = sort { $b->build_id <=> $a->build_id } $model->builds;
        for my $build (@builds) {
            last if $first_build and $second_build;
            my $build_revision = $build->software_revision;
            next unless defined $build_revision;
            if ($build_revision =~ /$first_revision/) {
                $first_build = $build;
            }
            elsif ($build_revision =~ /$second_revision/) {
                $second_build = $build;
            }
        }

        confess "Could not find build of model $model_id using revision $first_revision" unless $first_build;
        confess "Could not find build of model $model_id using revision $second_revision" unless $second_build;

        $self->status_message("Comparing build " . $first_build->build_id . " using revision $first_revision with build " . 
            $second_build->build_id . " using revision $second_revision from model $model_id");

        my %diffs = $first_build->compare_output($second_build->build_id);
        unless (%diffs) {
            $self->status_message("No differences found!");
        }
        else {
            my $diff_string = "DIFFERENCES FOUND:\n";
            for my $file (sort keys %diffs) {
                my $reason = $diffs{$file};
                $diff_string .= "  File: $file, Reason: $reason\n";
            }
            $self->status_message($diff_string);
        }
    }

    return 1;
}
1;  

