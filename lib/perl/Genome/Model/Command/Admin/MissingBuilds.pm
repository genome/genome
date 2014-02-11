package Genome::Model::Command::Admin::MissingBuilds;

class Genome::Model::Command::Admin::MissingBuilds {
    is => 'Genome::Command::Base',
    doc => 'Identify models that are missing builds.',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            is_optional => 1,
        },
        request_build => {
            is => 'Boolean',
            default => 0,
        },
    ],
};

use strict;
use warnings;
use Genome;

sub execute {
    my $self = shift;

    my @m = $self->models;
    unless (@m) {
        print STDERR "Using default search for models set to run as apipe-builder or ebelter...\n";
        @m = Genome::Model->get(
            'build_requested !=' => '1',
            run_as => ['apipe-builder', 'ebelter'],
        );
        print STDERR "Found " . @m . " models.\n";
    }

    my @m_need_build;
    for my $m (@m) {
        next if ($m->latest_build);
        next if ($m->build_requested && $m->build_requested =~ /^1$/);
        push @m_need_build, $m;
    }

    print STDERR "Found " . @m_need_build . " models that need a build.\n";

    for my $m (@m_need_build) {
        if ($self->request_build) {
            $m->build_requested(1);
        }
        else {
            print $m->id . "\t" . $m->processing_profile->name . "\n";
        }
    }

    if ($self->request_build) {
        print "Requested builds for " . @m_need_build . " models.\n";
    }
}
