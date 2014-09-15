package Genome::Test::Factory::AnalysisProject;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;

our @required_params = qw(name);

sub generate_obj {
    my $self = shift;
    my %params = @_;
    my $config_hash = delete $params{'config_hash'};

    unless(exists $params{status}) {
        $params{status} = 'In Progress'; #most commonly needed for tests
    }

    my $project = Genome::Config::AnalysisProject->create(%params);
    my $config_profile_item = Genome::Test::Factory::Config::Profile::Item->setup_object(analysis_project => $project);

    if ($config_hash) {
        for(values %$config_hash) {
            if(ref($_) eq 'ARRAY') {
                for(@$_){
                    $_->{config_profile_item} = $config_profile_item
                }
            } else {
                $_->{config_profile_item} = $config_profile_item
            }
        }
        $project->{__dummy_config_hash__} = $config_hash;
    }

    return $project;

}

sub create_name {
    return Genome::Test::Factory::Util::generate_name('analysis_project_name');
}

use Genome::Config::AnalysisProject;
my $old_get_config_method = \&Genome::Config::AnalysisProject::get_configuration_profile;

no warnings qw(redefine);
*Genome::Config::AnalysisProject::get_configuration_profile = sub {
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
    return UR::Util::deep_copy(shift->{'config'});
}

1;
