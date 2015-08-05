package Genome::Site::TGI::CleTest;

use strict;
use warnings;
use Genome;

class Genome::Site::TGI::CleTest {
    has_transient_optional => {
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
        },
    },
};
##NOTE: Do not change this file without approval from the CLE staff
sub get_builds {
    my $self = shift;
    my $cmd0 = Genome::Config::AnalysisProject::Command::Create->create(
        environment => "ad-hoc",
        name => "CLE Test",
        no_config => 1
    );
    my $analysis_project = $cmd0->execute;
    $self->analysis_project($analysis_project);
    my $menu_item1 = Genome::Config::AnalysisMenu::Item->get("9ab6e28f832a428393b87b171d444401");
    my $menu_item2 = Genome::Config::AnalysisMenu::Item->get("3770b8510d5a459f9c0bb01fabf56337");
    for my $tag_name (qw(discovery followup)) {
        my $tag = Genome::Config::Tag->get(name => $tag_name);
        Genome::Config::AnalysisProject::Command::AddMenuItem->execute(
            analysis_menu_items => $menu_item1,
            tags => $tag,
            analysis_project => $analysis_project,
        );
    }
    for my $tag_name (qw(germline)) {
        my $tag = Genome::Config::Tag->get(name => $tag_name);
        Genome::Config::AnalysisProject::Command::AddMenuItem->execute(
            analysis_menu_items => $menu_item2,
            tags => $tag,
            analysis_project => $analysis_project,
        );
    }
    my $subject_mapping_file = Genome::Sys->create_temp_file_path;
    my $subject_mapping = <<'SUBJECT_MAPPINGS';
H_KA-174556-1309237	H_KA-174556-1309246				followup
H_KA-174556-1309245	H_KA-174556-1309246				discovery
H_KA-174556-1309246					germline
SUBJECT_MAPPINGS

    Genome::Sys->write_file($subject_mapping_file, $subject_mapping);
    Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation->execute(
        analysis_project => $analysis_project,
        file_path => $subject_mapping_file,
    );

    my @instrument_data = qw(2893814999 2893815000 2893815001 2893815002 2893815003 2893815016 2893815018 2893815020 2893815023 2893815024 2893815447 2893815448 2893815449 2893815451 2893815484 2893815485 2893815489 2893815492);
    for my $id (@instrument_data) {
        my $instrument_data = Genome::InstrumentData->get($id);
        Genome::Config::AnalysisProject::InstrumentDataBridge->create(
            analysis_project => $analysis_project,
            instrument_data => $instrument_data,
            status => 'new',
        );
    }

    Genome::Config::AnalysisProject::Command::Release->execute(
        analysis_projects => [$analysis_project],
    );

    my @instrument_data_objects = Genome::InstrumentData->get(id => \@instrument_data);
    Genome::Config::Command::ConfigureQueuedInstrumentData->execute(
        instrument_data => \@instrument_data_objects,
    );

    my @models = Genome::Model->get(analysis_project => $analysis_project);
    Genome::Model::Build::Command::Start->execute(
        models => \@models,
    );
    UR::Context->commit;
    my @builds = Genome::Model::Build->get(model_id => [map {$_->id} @models]);
    return @builds;
}

sub get_process {
    my $self = shift;
    my $analysis_project = $self->analysis_project;
    unless (defined $analysis_project) {
        die $self->error_message("Must run get_builds first");
    }

    my $discovery_subject = $self->get_sample_from_subject_mapping('discovery');
    my $followup_subject = $self->get_sample_from_subject_mapping('followup');
    my $germline_subject = $self->get_sample_from_subject_mapping('germline');

    my @models_for_process = Genome::Model::SomaticValidation->get(analysis_project => $analysis_project,
        region_of_interest_set_name => 'SeqCap EZ Human Exome v3.0 + AML RMG pooled probes + WO2830729 pooled probes + WO2840081 pooled probes');
    my @coverage_models_for_process = Genome::Model::SomaticValidation->get(analysis_project => $analysis_project,
        tumor_sample => [$discovery_subject, $followup_subject]);
    my @germline_models = grep {$_->tumor_sample eq $germline_subject} @models_for_process;
    Genome::VariantReporting::Command::Wrappers::Trio->execute(
        models => \@models_for_process,
        coverage_models => \@coverage_models_for_process,
        tumor_sample => $discovery_subject,
        followup_sample => $followup_subject,
        normal_sample => $germline_subject,
    );
    UR::Context->commit;
    return $process;
}

sub get_sample_from_subject_mapping {
    my ($self, $name) = @_;
    my $tag = Genome::Config::Tag->get(name => $name);
    my $mapping = Genome::Config::AnalysisProject::SubjectMapping->get(analysis_project => $self->analysis_project,
        tags => [$tag]);
    my ($bridge) = $mapping->subject_bridges(label => 'tumor_sample');
    return $bridge->subject;
}

sub diff_build {
    my $self = shift;
    my $build = shift;

    printf("Starting diff (new build = %s)...\n", $build->id);

    my @blessed_builds = Genome::Model::Build->get(id => [qw(
    185d8bac3d7c4437b7ce9207dfadac7e
    5f57f08a8fd84886b275f0b5571e9fd7
    e02e8a5ccaad458d839de51eb47d8d2c
    26e65adaa8034dd99ef92b27f61ad862
    3243f261a8b64c089a5254291f7c2de3
    f87701c292e843958d088e171a65a67a
    301d5d51c96e41308d015008720f3962
    )]);

    my $matching_blessed_build;
    for my $blessed_build (@blessed_builds) {
        if ($build->tumor_sample eq $blessed_build->tumor_sample and
                ((!defined($build->normal_sample) and !defined($blessed_build->normal_sample)) or
                    ($build->normal_sample eq $blessed_build->normal_sample)) and
            $build->target_region_set_name eq $blessed_build->target_region_set_name) {
                $matching_blessed_build = $blessed_build;
                last;
            }
    }
    unless (defined $matching_blessed_build) {
        die $self->error_message("No matching blessed build found for build %s\n", $build->id);
    }

    my $diff_cmd = Genome::Model::Build::Command::Diff->create(
        new_build => $build,
        blessed_build => $matching_blessed_build,
    );
    unless ($diff_cmd->execute) {
        die $self->error_message("Diff command failed to execute for build %s!\n", $build->id);
    }

    if ($diff_cmd->has_diffs) {
        return $diff_cmd;
    }
    return;
}

sub diff_process {
    my $self = shift;
    my $process = shift;

    my $blessed_process = Genome::Process->get("2b2b3d8481284fcf8a1632d65ba58083");

    printf("Starting diff (new process = %s)...\n", $process->id);
    my $diff_cmd = Genome::Process::Command::Diff->create(
        new_process => $process,
        blessed_process => $blessed_process,
    );
    unless ($diff_cmd->execute) {
        die $self->error_message("Diff command failed to execute!\n");
    }

    if ($diff_cmd->has_diffs) {
        return $diff_cmd;
    }
    return;
}
1;
