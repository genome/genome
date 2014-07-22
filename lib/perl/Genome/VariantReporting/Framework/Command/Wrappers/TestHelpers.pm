package Genome::VariantReporting::Framework::Command::Wrappers::TestHelpers;

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::Model::SomaticValidation;
use Genome::Test::Factory::ProcessingProfile::SomaticValidation;
use Genome::Test::Factory::Build;
use Exporter 'import';

our @EXPORT_OK = qw(get_build);

my $pp = Genome::Test::Factory::ProcessingProfile::SomaticValidation->setup_object();

my $start_id = "-b8bc947b2fbf464592508cc021fa63ef";
my $fl_counter;
sub _get_or_create_feature_list {
    my $name = shift;
    my $feature_list = Genome::FeatureList->get(name => $name);
    unless ($feature_list) {
        my $roi = Genome::FeatureList->__define__(name => $name, id => $start_id.$fl_counter);
        $fl_counter++;
    }
    return $feature_list;
}

sub get_build {
    my ($roi_name, $tumor_sample, $normal_sample) = @_;
    _get_or_create_feature_list($roi_name);
    my $discovery_model = Genome::Test::Factory::Model::SomaticValidation->setup_object(processing_profile_id => $pp->id);
    $discovery_model->tumor_sample($tumor_sample);
    $discovery_model->normal_sample($normal_sample);
    my $discovery_build = Genome::Test::Factory::Build->setup_object(model_id => $discovery_model->id);
    $discovery_model->region_of_interest_set_name($roi_name);
    $discovery_build->region_of_interest_set_name($roi_name);
    #$discovery_build->tumor_sample($tumor_sample);
    #$discovery_build->normal_sample($normal_sample);
    return $discovery_build;
}

1;

