#!/usr/bin/perl

# This script grabs various objects needed for the cle-test
# If this script changes, you will need to create a new template for
# the Cle Regression Test jenkins job:
# 1. Use genome-env -d template1 to drop into an environment with a clean testdb
# 2. Set up the test_filldb config to point to the clean db (see Genome::DataSource::GMSchema)
# 3. Set up the gmschema configs to point to the production db (see Genome::DataSource::GMSchema)
# 4. Run this script
# 5. use `test-db template create` to create a new template
# 6. Use the new template in the jenkins configuration
use Genome;
use Genome::Site::TGI::CleTest;
use List::MoreUtils qw(uniq);

#things we need to load
# Instrument data and its attributes
# Diff models
# Annotation
# Reference sequence
# Reference sequence converters
# sample attributes
# feature lists
# variant report to diff
# apipe tester user object
# analysis project menu items
# Genome::Config::Tag items.

sub load_instrument_data {
    my @instrument_data_ids = @_;

    my @inst_data = Genome::InstrumentData->get(id => \@instrument_data_ids);
    for my $inst_data_obj (@inst_data) {
        load_through_methods($inst_data_obj, qw( library attributes ));
        my $lib = $inst_data_obj->library; #libraries have no attributes
        load_through_methods($inst_data_obj->individual, 'attributes');
        load_through_methods($inst_data_obj->sample, 'attributes');
        load_through_methods($inst_data_obj->taxon, 'attributes');
    }
}

sub load_builds {
    my @build_ids = @_;
    load_builds_without_results(@build_ids);
    my @builds = Genome::Model::Build->get(id => \@build_ids);
    for my $build (@builds) {
        load_through_methods($build, qw( results ));
    }
}

sub load_imported_annotation_builds {
    my @build_ids = @_;
    load_builds(@build_ids);
    my @builds = Genome::Model::Build->get(id => \@build_ids);
    for my $build (@builds) {
        load_through_methods($build, qw( snv_result indel_result));
    }
}

sub load_processing_profiles {
    my @processing_profile_ids = @_;
    my @processing_profiles = Genome::ProcessingProfile->get(\@processing_profile_ids);
    for my $pp (@processing_profiles) {
        load_through_methods($pp, qw( params ));
    }
}

sub load_processes {
    my @process_ids = @_;
    my @processes = Genome::Process->get(\@process_ids);
    for my $process (@processes) {
        load_through_methods($process, qw( disk_allocation results inputs ));
    }
}

sub load_builds_without_results {
    my @build_ids = @_;
    my @builds = Genome::Model::Build->get(id => \@build_ids);
    for my $build (@builds) {
        load_through_methods($build, qw(inputs metrics));
        load_through_methods($build->model, qw( subject inputs ));
        load_through_methods($build->model->processing_profile, qw( params ));
    }
}

sub load_software_results {
    my @software_results = @_;
    for my $result (@software_results) {
        my @inputs = $result->inputs;
        my @params = $result->params;
        my @metrics = $result->metrics;
        my @disk_allocations = $result->disk_allocations;
    }
}

sub load_through_methods {
    my ($obj, @methods) = @_;
    for my $method (@methods) {
        my @results = $obj->$method; #should probably error handle this crap somehow
        for my $result (@results) {
            if(defined $result and $result->isa('Genome::SoftwareResult')) {
                load_software_results($result);
            }
        }
    }
}
    
sub load_converters {
    my @reference_ids = @_;
    for my $from (@reference_ids) {
        for my $to (@reference_ids) {
            if($from != $to) {
                my @converters = Genome::Model::Build::ReferenceSequence::Converter->get( source_reference_build_id => $from, destination_reference_build_id => $to );
                for my $converter (@converters) {
                #die "More than one converter found from: $from to: $to!\n" if (@converters > 1);

                    load_software_results($converter) if $converter;
                }
            }
        }
    }
}

sub load_feature_lists {
    my @feature_list_ids = @_;
    
    my @feature_lists = Genome::FeatureList->get(\@feature_list_ids);
    for my $list (@feature_lists) {
        load_through_methods($list, qw( disk_allocation ));
#        load_builds($list->reference_id); #ASSUMING YOU HAVE DONE THE WORK TO KNOW YOUR REFERENCE SEQUENCE
    }
}

sub load_menu_items {
    my @ids = @_;
    my @menu_items = Genome::Config::AnalysisMenu::Item->get(\@ids);
}

sub load_config_tags_by_name {
    my @names = @_;
    my @tags = Genome::Config::Tag->get( name => \@names );
}

sub load_users_by_username {
    my @usernames = @_;
    my @users = Genome::Sys::User->get(username => \@usernames);
}

sub main {

    my $config = Genome::Site::TGI::CleTest::get_config();
    my @instrument_data_ids = @{$config->{instrument_data}};
    load_instrument_data(@instrument_data_ids);

    #$ genome model list analysis_projects.id=b0db77fc61334bc69527a7000b643ee6,subject.name~H_KA-174556% -sh last_succeeded_build.id --noheader
    my @blessed_builds = @{$config->{blessed_builds}};
    load_builds(@blessed_builds);

    my @build_ids = @{$config->{annotation_builds}};

    my @imported_variation_list_ids = @{$config->{imported_variation_builds}};

    load_builds(@build_ids);
    load_imported_annotation_builds(@imported_variation_list_ids);

    my @reference_sequences = @{$config->{reference_sequence_builds}};
    load_builds_without_results(@reference_sequences);
    load_converters(@reference_sequences);


    my @featurelist_ids = @{$config->{feature_lists}};
    load_feature_lists(@featurelist_ids);

    my %tag_to_menu_item = %{$config->{tag_to_menu_item}};
    my @tag_names = keys(%tag_to_menu_item);
    load_config_tags_by_name(@tag_names);
    my @menu_items = uniq(values %tag_to_menu_item);
    load_menu_items(@menu_items);

    my @processing_profiles = @{$config->{processing_profiles}};
    load_processing_profiles(@processing_profiles);

    load_users_by_username('apipe-tester');

    my $process = $config->{blessed_process};
    load_processes($process); #our report process to diff

    my @misc_software_results = Genome::SoftwareResult->get(
        id => $config->{misc_software_results});
    load_software_results(@misc_software_results);
    return UR::Context->commit();
}

main();

