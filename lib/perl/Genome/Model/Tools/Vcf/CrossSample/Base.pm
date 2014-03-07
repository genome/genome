package Genome::Model::Tools::Vcf::CrossSample::Base;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::Sam;

class Genome::Model::Tools::Vcf::CrossSample::Base {
    is => ['Command::V2', 'Genome::Model::Tools::Park::Base'],
    is_abstract => 1,

    has_input => [
        output_directory => {
            is => 'Text',
            is_optional => 0,
            doc => 'the directory where you want results stored',
        },
    ],
    has_optional_input => [
        builds => {
            is => 'Genome::Model::Build',
            require_user_verify => 0,
            is_many => 1,
            doc => 'The builds that you wish to create a cross-sample vcf for',
        },
        model_group => {
            is => 'Genome::ModelGroup',
            require_user_verify => 0,
            doc => 'Model group from which last succeeded builds will be pulled',
        },
        allow_multiple_processing_profiles => {
            is => 'Boolean',
            default => 0,
            doc => 'Setting this prevents the check for identical processing profiles on all input builds',
        },

        # region limiting
        roi_list => {
            is => 'Genome::FeatureList',
            doc => 'Set this to limit the incoming vcfs to roi target regions',
        },
        wingspan => {
            is => 'Text',
            doc => 'Set this to add a wingspan to region limiting',
        },

        joinx_version => {
            is => 'Text',
            doc => "Version of joinx to use, will be resolved to the latest default if not specified",
        },
        merge_strategy_file => {
            is => 'Path',
            default => 'NONE',
            doc => "The merge strategy file, provided to joinx to specify how to handle info field merging (-M option)",
        },

        samtools_version => {
            is => 'Text',
            doc => "Version of samtools to use, will be resolved to the latest default if not specified",
        },
        samtools_minimum_base_quality => {
            is => 'Number',
            default => 13,
            doc => "Minimum base quality for a base to be considered.",
        },
        samtools_minimum_mapping_quality => {
            is => 'Number',
            default => 10,
            doc => "Minimum mapping quality for an alignment to be used.",
        },
    ],
    has_optional_transient => [
        _reference_sequence_build => {
            is => 'Genome::Model::Build',
        },
        reference => {
            is => 'Path',
        },
        reference_index => {
            is => 'Path',
        },
        vcfs => {
            is => 'Path',
            is_many => 1,
        },
    ],
};


sub help_detail {
    return <<EOS
We first get the appropriate vcf file from the last_succeeded_build in
each model in the model-group.  From those we create a vcf that includes
every variation found (in any vcf).  For every location in this vcf we look
to the bam file for evidence (this is done for every build separately).
These 'back-filled' vcfs are then combined into a single multi-sample vcf.
EOS
}

sub vcf_accessor {
    die "Abstract";
}

sub _execute {
    # needs to return the process uri
    die "Abstract";
}

sub execute {
    my $self = shift;

    $self->_resolve_joinx_version;
    $self->_resolve_samtools_version;
    $self->_resolve_builds;
    $self->_validate_inputs;

    $self->status_message("Running Workflow...");
    my $process_uri = $self->_execute;

    $self->status_message("Symlinking Results...");
    $self->_symlink_results($process_uri);

    return 1;
}

sub _resolve_joinx_version {
    my $self = shift;
    unless (defined $self->joinx_version) {
        $self->joinx_version(Genome::Model::Tools::Joinx->get_default_version);
    }
    unless ($self->joinx_version >= 1.8) {
        die $self->error_message("Joinx version 1.8 or greater is required");
    }
    return;
}

sub _resolve_samtools_version {
    my $self = shift;
    unless (defined $self->samtools_version) {
        $self->samtools_version(Genome::Model::Tools::Sam::default_samtools_version);
    }
    return;
}

sub _resolve_builds {
    my $self = shift;

    $self->status_message("Resolving Builds...");
    my @builds;
    if ($self->builds and not $self->model_group) {
        @builds = $self->builds;
    } elsif ($self->model_group and not $self->builds) {
        my $command = Genome::ModelGroup::Command::GetLastCompletedBuilds->execute(
            model_group => $self->model_group);
        @builds = $command->builds;
        $self->builds(\@builds);
    }
    else {
        die $self->error_message("Given both builds and model-groups or neither.");
    }

    return;
}

sub _validate_inputs {
    my $self = shift;

    $self->status_message("Validating Inputs...");

    $self->_validate_builds;
    $self->_check_build_files;
    return;
}

sub _roi_bed_file {
    my $self = shift;
    return $self->roi_list->resolve_bed_for_reference($self->_reference_sequence_build);
}

sub _roi_name {
    my $self = shift;
    return $self->roi_list->name;
}

sub _validate_builds {
    my $self = shift;

    my @builds = $self->builds;

    my $first_build = $builds[0];
    my $ref = $first_build->reference_sequence_build;
    $self->_reference_sequence_build($ref);
    $self->reference($ref->full_consensus_path('fa'));
    $self->reference_index($ref->full_consensus_path('fa.fai'));
    my $pp = $first_build->processing_profile;
    my %validation_params = (
        builds => \@builds,
        builds_can => [qw(reference_sequence_build whole_rmdup_bam_file get_snvs_vcf get_indels_vcf)],
        status => ['Succeeded'],
        reference_sequence => [$ref],
    );
    if (!$self->allow_multiple_processing_profiles) {
        $validation_params{'processing_profile'} = [$pp];
    }
    Genome::Model::Build::Command::Validate->execute(%validation_params);
    return;
}

sub _check_build_files {
    my $self = shift;

    my (@builds_with_file, @builds_without_file, @vcf_files);
    my $accessor = $self->vcf_accessor;
    for my $build ($self->builds) {
        my $vcf_file = $build->$accessor;
        if (-s $vcf_file) {
            push @builds_with_file, $build->id;
            push @vcf_files, $vcf_file;
        } else {
            push @builds_without_file, $build->id;
        }
    }
    $self->vcfs(\@vcf_files);

    my $num_builds = scalar($self->builds);
    unless( scalar(@builds_with_file) == $num_builds){
        die $self->error_message("The number of input builds ($num_builds) did not match the" .
            " number of vcf files found (" . scalar (@builds_with_file) . ").\n" .
            "Check the input builds for completeness.\n" .
            "Builds with a file present: " . join(",", @builds_with_file) . "\n" .
            "Builds with missing or zero size file: " . join(",", @builds_without_file) . "\n"
        );
    }
    return;
}

sub _symlink_results {
    my ($self, $process_uri) = @_;

    $self->_link_process($process_uri, $self->output_directory);
    return;
}

sub sample_names {
    my $self = shift;
    my @samples;
    for my $build ($self->builds) {
        push @samples, $build->subject->name;
    }
    return @samples;
}

1;
