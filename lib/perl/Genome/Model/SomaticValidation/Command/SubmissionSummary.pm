package Genome::Model::SomaticValidation::Command::SubmissionSummary;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::SomaticValidation::Command::SubmissionSummary {
    is => 'Genome::Command::Base',
    doc => "List a summary of the merged alignment BAMs for the provided builds and a file suitable for submitting the bam list.",
    has => [
        sample_mapping_file => {
            is=> 'String',
            doc => 'this command will generate this list of bam file names & samples',
        }, 
        bam_list_file => {
            is=> 'String',
            doc => 'this command will generate this list of bam file names, newline separated suitable for gxfer to submit bams',
        },
        md5_list_file => {
            is=> 'String',
            doc => 'this will generate a list of bam md5 file names, newline separated',
        }
    ],
    has_optional => [
        builds => {
            is => 'Genome::Model::Build::SomaticValidation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'List of builds to use if searching on a list of builds'
        },
        flow_cell_id => {
            is => 'String',
            shell_args_position => 2,
            doc => 'Flow Cell ID to use if searching on a flowcell'
        },
        flow_cell_subset_name => {
            is => 'String',
            shell_args_position => 3,
            doc=> 'Flow cell subset name to use if searching on a flow cell'
        },
        reference_sequence_name => {
            is=>'String',
            doc=>'Only include models with this reference sequence name',
        },
        region_of_interest_set_name => {
            is=>'String',
            doc=>'Only include models with this region of interest set name',
        },
        exclude => {
            is=>'String',
            doc=>'Don\'t include models that contain this string in the name (ie. "Pooled_Library")',
        }
    ],
};


sub help_detail {
    return "List the path of the merged alignment BAMs for the provided builds/flow cells, and generate a sample mapping.";
}


sub execute {
    my $self = shift;

    die "SomaticValidation models don't create md5 checksums. This tool cannot yet be run.\n";

    if ($self->flow_cell_id && $self->builds) {
        $self->error_message("Ambiguous input; you must provide either a flow cell id or a set of builds -- not both. ");
        return;
    }

    if ($self->flow_cell_id) {
        my @instr_data_params = (flow_cell_id=>$self->flow_cell_id);
        push @instr_data_params, (subset_name => $self->flow_cell_subset_name) if $self->flow_cell_subset_name;
        $self->status_message(sprintf("Searching for instrument data for %s/%s...  ", $self->flow_cell_id, defined $self->flow_cell_subset_name ? $self->flow_cell_subset_name : "all" ));
        my @instr_data = Genome::InstrumentData::Solexa->get(@instr_data_params);
        $self->status_message(sprintf("Found %s instrument data.\n", scalar @instr_data));

        my %model_ids = map {$_->model_id, 1} map {Genome::Model::Input->get(name=>'instrument_data', value_id=>$_->id)} @instr_data;

        my @raw_models = Genome::Model->get(id=>[keys %model_ids]);
        my @models;

        # filter out Lane QC models, Pooled_Library models, and only the ROI/refseq requested if there was one
        for (@raw_models) {
            push @models, $_ unless (($_->subject_name =~ m/^Pooled_Library/) ||
                                     ($_->processing_profile->append_event_steps && $_->processing_profile->append_event_steps =~ m/LaneQc/) ||
                                     ($self->region_of_interest_set_name && $_->region_of_interest_set_name && $_->region_of_interest_set_name ne $self->region_of_interest_set_name) ||
                                     ($self->reference_sequence_name && $_->reference_sequence_build->name ne $self->reference_sequence_name));

        }

        $self->status_message(sprintf("Found %s models", scalar @models));

        my @builds;
        for (@models) {
            my ($latest_build) = sort {$b->id <=> $a->id} grep {$_->status eq 'Succeeded'} $_->builds;
            push @builds, $latest_build if $latest_build;
        }

        $self->builds([@builds]);

    } elsif ($self->builds) {
        $self->status_message("Using user-supplied set of builds...");

    } else {
        $self->error_message("You must provide either a flow cell id or a set of builds.  ");
        return;
    }

    my @filtered_builds;
    for my $b ($self->builds) {
        if ($b->status ne "Succeeded") {
            warn sprintf("Filtering out build %s (model %s) because its status is %s, not succeeded.", $b->id, $b->model->name, $b->status);
        } else {
            push @filtered_builds, $b;
        }
    }
    $self->builds([@filtered_builds]);

    my $samp_map = IO::File->new(">".$self->sample_mapping_file);
    unless ($samp_map) {
        $self->error_message("Failed to open sample mapping file for writing ". $self->sample_mapping_file);
        return;
    }
    my $bam_list = IO::File->new(">".$self->bam_list_file);
    unless ($bam_list) {
        $self->error_message("Failed to open bam list file for writing ". $self->bam_list_file);
        return;
    }

    my $md5_list = IO::File->new(">".$self->md5_list_file);
    unless ($md5_list) {
        $self->error_message("Failed to open md5 list file for writing ". $self->md5_list_file);
        return;
    }


    my @builds = grep {$_ == $_->model->last_succeeded_build} $self->builds;

    # If we have models to exclude by name, do so now...in 5...4...3...2...1
    my $exclude = $self->exclude;
    @builds = grep {not $_->model->name =~ /$exclude/} @builds if $exclude;

    for my $build (@builds) {
        my $roi_name = $build->model->region_of_interest_set_name ? $build->model->region_of_interest_set_name : 'N/A';
        my $refbuild_name = $build->model->reference_sequence_build->name ? $build->model->reference_sequence_build->name : 'N/A';

        print $samp_map join ("\t", $build->tumor_sample->name, $refbuild_name, $roi_name, $build->tumor_bam, basename($build->tumor_bam));
        print $samp_map "\n";
        if($build->normal_sample) {
            print $samp_map join ("\t", $build->normal_sample->name, $refbuild_name, $roi_name, $build->normal_bam, basename($build->normal_bam));
            print $samp_map "\n";
        }
    }

    my @bams = map {$_->tumor_bam} @builds;
    push @bams, grep { defined $_ } map {$_->normal_bam} @builds;
    print $bam_list join ("\n", @bams);
    print $md5_list join ("\n", map {$_.".md5"} @bams);

    $samp_map->close;
    $bam_list->close;
    $md5_list->close;

    return 1;
}

1;

