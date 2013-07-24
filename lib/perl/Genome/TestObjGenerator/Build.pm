package Genome::TestObjGenerator::Build;
use Genome::TestObjGenerator::Base;
@ISA = (Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;

my @required_params = ("model_id", "data_directory");

sub generate_obj {
    my $self = shift;
    my %params = @_;
    my $status = delete $params{status};
    my $build = Genome::Model::Build->create(%params);
    if (defined $status) {
        _set_status_on_build($build, $status);
    }
    return $build;
}

sub get_required_params {
    return \@required_params;
}

sub create_data_directory {
    return Genome::Sys->create_temp_directory();
}

sub _set_status_on_build {
    my $build = shift;
    my $status = shift;

    $build->the_master_event->event_status("Succeeded");
}

1;
