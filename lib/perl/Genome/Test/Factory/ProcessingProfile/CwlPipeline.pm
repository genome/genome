package Genome::Test::Factory::ProcessingProfile::CwlPipeline;
use Genome::Test::Factory::ProcessingProfile;
@ISA = (Genome::Test::Factory::ProcessingProfile);

our @required_params = qw(name main_workflow_file primary_docker_image);

sub create_main_workflow_file {
    return '/dev/null/main.yml';
}

sub create_primary_docker_image {
    return 'mgibio/docker-test';
}

1;
