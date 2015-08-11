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
                die "More than one converter found from: $from to: $to!\n" if (@converters > 1);

                load_software_results($converters[0]) if @converters;
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
    my @instrument_data_ids = qw(
        2893814999
        2893815000
        2893815001
        2893815002
        2893815003
        2893815016
        2893815018
        2893815020
        2893815023
        2893815024
        2893815447
        2893815448
        2893815449
        2893815451
        2893815484
        2893815485
        2893815489
        2893815492
    );
    load_instrument_data(@instrument_data_ids);

    #$ genome model list analysis_projects.id=b0db77fc61334bc69527a7000b643ee6,subject.name~H_KA-174556% -sh last_succeeded_build.id --noheader
    my @blessed_builds = qw(
    185d8bac3d7c4437b7ce9207dfadac7e
    5f57f08a8fd84886b275f0b5571e9fd7
    e02e8a5ccaad458d839de51eb47d8d2c
    26e65adaa8034dd99ef92b27f61ad862
    3243f261a8b64c089a5254291f7c2de3
    f87701c292e843958d088e171a65a67a
    301d5d51c96e41308d015008720f3962
    );
    load_builds(@blessed_builds);

    my @build_ids = (
        'd00a39c84382427fa0efdec3229e8f5f', #annotation build
    );

    my @imported_variation_list_ids = (
        '7fd9ca9a6bf74a098c23b32c1b00fbc0', #dbsnp
        '126747180', #NHLBI version 2012.07.23
        '847b3cacad1249b8b7e46f89e02d96da', #DOCM SNV
        '1040bf09070c4176a5256fa8a075378f', #DOCM INDEL

    );
    load_builds(@build_ids);
    load_imported_annotation_builds(@imported_variation_list_ids);
    
    my @reference_sequences = (
        '106942997', #build37, other references may crop up in feature lists
        '108563338', #nimblegen refseq for feature list
        '102671028', #g1k-human-build37 for segdups
    );
    load_builds_without_results(@reference_sequences);
    load_converters(@reference_sequences);


    my @featurelist_ids = (
        'a2733b8111c24110805025137c878e13', #Full ROI
        '0e4973c600244c3f804d54bee6f81145', #RMG ROI
        'ef6583623a614c8ebaea471a4c9748fc', #AML Complex Mutation Region ROI
        '696318bab30d47d49fab9afa845691b', #homopolymer regions from G::VR::C::Wrapper::ModelPair
        '4e77d10f29f44a0792ebc8d6ea0c4a2b', #segdups my $b = Genome::Model::Build->get('26e65adaa8034dd99ef92b27f61ad862'); $b->get_feature_list("segmental_duplications")->id

    );
    load_feature_lists(@featurelist_ids);

    my @tag_names = qw( discovery followup germline );
    load_config_tags_by_name(@tag_names);
    my @menu_items = (
        '3770b8510d5a459f9c0bb01fabf56337',#germline
        '9ab6e28f832a428393b87b171d444401',#somatic
    );
    load_menu_items(@menu_items);

    my @processing_profiles = (
        'a813f8067a8c4c9791fec53dd29d85ca', #somatic CLE PP
        '6cab54acdf704c7c9e8adc7aa8facff4', #germline CLE PP
    );
    load_processing_profiles(@processing_profiles);

    load_users_by_username('apipe-tester');

    load_processes('2b2b3d8481284fcf8a1632d65ba58083'); #our report process to diff

    return UR::Context->commit();
}

main();

