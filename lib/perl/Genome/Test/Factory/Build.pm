package Genome::Test::Factory::Build;
use Genome::Test::Factory::Base;
@ISA = (Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;

our @required_params = qw(model_id data_directory run_by status);

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

sub create_data_directory {
    return Genome::Sys->create_temp_directory();
}

sub _set_status_on_build {
    my $build = shift;
    my $status = shift;

    $build->status($status);
    $build->date_completed("2013-07-11 20:47:51");
}

sub create_run_by {
    return Genome::Sys->username;
}

sub create_status {
    return 'New';
}

1;
