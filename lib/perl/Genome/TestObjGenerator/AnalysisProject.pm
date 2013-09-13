package Genome::TestObjGenerator::AnalysisProject;
use base qw(Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;

our @required_params = qw(name);

sub generate_obj {
    my $self = shift;
    my %params = @_;
    my $config_hash = delete $params{'config_hash'};
    my $project = Genome::Config::AnalysisProject->create(%params);

    if ($config_hash) {
        $project->{__dummy_config_hash__} = $config_hash;
    }

    return $project;

}

sub create_name {
    return Genome::TestObjGenerator::Util::generate_name('analysis_project_name');
}

use Genome::Config::AnalysisProject;
my $old_get_config_method = \&Genome::Config::AnalysisProject::get_configuration_reader;

no warnings qw(redefine);
*Genome::Config::AnalysisProject::get_configuration_reader = sub {
    my $self = shift;
    if ($self->{__dummy_config_hash__}) {
        return bless { config => $self->{__dummy_config_hash__} },
            'Genome::Config::DummyConfigReader';
    } else {
        return $old_get_config_method->($self, @_);
    }
};
use warnings;

package Genome::Config::DummyConfigReader;
sub get_config {
    return shift->{'config'};
}

1;
