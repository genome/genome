package Genome::Model::CwlPipeline::Runner::Toil;

use strict;
use warnings;

use File::Spec;

use Genome;

class Genome::Model::CwlPipeline::Runner::Toil {
    is => 'UR::Singleton',
};

sub version_exists {
    my $class = shift;
    my $version = shift;

    my $allocation = Genome::Disk::Allocation->get(
        allocation_path => File::Spec->join("toil", $version),
        status => 'active',
    );

    return !! $allocation;
}

sub cwl_runner_wrapper_for_version {
    my $class = shift;
    my $version = shift;

    my $allocation = Genome::Disk::Allocation->get(
        allocation_path => File::Spec->join("toil", $version),
        status => 'active',
    );
    unless ($allocation) {
        $class->fatal_message('No toil version found for version %s', $version);
    }

    my $wrapper = File::Spec->join($allocation->absolute_path, 'toil-cwl-runner-wrapper.sh');
    unless (-e $wrapper) {
        $class->fatal_message('Wrapper script not found in %s for version %s', $allocation->absolute_path, $version);
    }

    return $wrapper;
}

1;
