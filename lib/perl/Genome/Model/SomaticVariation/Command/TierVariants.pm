package Genome::Model::SomaticVariation::Command::TierVariants;

use strict;
use warnings;
use Genome;

class Genome::Model::SomaticVariation::Command::TierVariants{
    is => 'Genome::Command::Base',
    has =>[
        build_id => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticVariation model',
        },
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
            id_by => 'build_id',
        },
        _tier_file_location => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
    ],
};

sub shortcut {
    my $self = shift;
    my $build = $self->build;

    my @existing_results;
    for my $qual ('hq', 'lq') {
        TYPE: for my $variant_type ('snv', 'indel') {
            my $accessor = $variant_type . '_detection_strategy';
            next TYPE unless $build->$accessor;

            my @params = $self->params_for_result($qual, $variant_type);
            return unless @params; #can't generate result so cannot shortcut

            my $existing_result = Genome::Model::Tools::DetectVariants2::Classify::Tier->get_with_lock(@params);
            if($existing_result) {
                push @existing_results, $existing_result;
            } else {
                return;
            }
        }
    }

    for my $existing_result (@existing_results) {
        $self->link_result_to_build($existing_result);
    }
}


sub execute {
    my $self = shift;
    my $build = $self->build;
    unless ($build){
        die $self->error_message("no build provided!");
    }

    $self->debug_message("executing tier variants step on snvs and indels");

    for my $qual ('hq', 'lq') {
        TYPE: for my $variant_type ('snv', 'indel') {
            my $accessor = $variant_type . '_detection_strategy';
            next TYPE unless $build->$accessor;

            my @params = $self->params_for_result($qual, $variant_type);
            if(@params) {
                my $result = Genome::Model::Tools::DetectVariants2::Classify::Tier->get_or_create(@params);
                if($result) {
                    $self->link_result_to_build($result);
                } else {
                    die $self->error_message('Failed to generate result for ' . $qual . ' ' . $variant_type);
                }
            } else {
                $self->_tier_non_result_files($qual, $variant_type);
            }
        }
    }

    $self->debug_message("Tier Variants step completed");
    return 1;
}

sub _tier_non_result_files {
    my $self = shift;
    my $qual = shift;
    my $variant_type = shift;
    my $build = $self->build;

    #TODO find a better way to set this. It should perhaps come from the build?
    my $bed_version = 2;

    my $tiering_version = $build->tiering_version;
    $self->debug_message("Using tiering_bed_files version ".$tiering_version);
    my $tier_file_location = $build->annotation_build->tiering_bed_files_by_version($tiering_version);

    unless (-d $tier_file_location){
        die $self->error_message("Couldn't find tiering bed files from annotation build");
    }

    $self->_tier_file_location($tier_file_location);

    my $effects_dir = $self->build->data_directory."/effects";
    unless(-d $effects_dir){
        Genome::Sys->create_directory($effects_dir);
    }

    my @name_sets;
    if($qual eq 'hq') {
        @name_sets = (['novel', $variant_type . 's.hq.novel'], ['novel', $variant_type . 's.hq.previously_detected']);
    } elsif ($qual eq 'lq') {
        @name_sets = (['variants', $variant_type . 's.lq']);
    } else {
        die ('unexpected qual ' . $qual);
    }

    for my $name_set (@name_sets){
            $self->run_fast_tier($name_set,$bed_version, $tiering_version, 'bed');
    }

    return 1;
}


sub run_fast_tier {
    my $self = shift;
    my ($name_set, $bed_version, $tiering_version, $format) =  @_;
    #breaking up filename and subdir parts of the data set path so we can put the output tiering files in the effects directory
    my $dir = $$name_set[0];
    my $name = $$name_set[1];

    my $build = $self->build;

    my $path_to_tier = $build->data_set_path("$dir/$name",$bed_version,$format);
    unless (-e $path_to_tier){
        die $self->error_message("No $name file for build!");
    }

    my ($tier1_path, $tier2_path, $tier3_path, $tier4_path) = map {$build->data_set_path("effects/$name.tier".$_,$bed_version,$format)}(1..4);

    my %params;

    #Skip line count on fast-tiering if running on input with duplicates (lq, in this case)
    my $lq = $name =~ m/lq/;

    if (-s $path_to_tier){
        %params = (
            variant_bed_file => $path_to_tier,
            tier_file_location => $self->_tier_file_location,
            tiering_version => $tiering_version,
            skip_line_count => $lq,
            tier1_output => $tier1_path,
            tier2_output => $tier2_path,
            tier3_output => $tier3_path,
            tier4_output => $tier4_path,
        );
        my $tier_snvs_command = Genome::Model::Tools::FastTier::FastTier->create(%params);
        unless ($tier_snvs_command){
            die $self->error_message("Couldn't create fast tier command from params:\n" . Data::Dumper::Dumper(\%params));
        }
        my $snv_rv = $tier_snvs_command->execute;
        my $snv_err =$@;
        unless($snv_rv){
            die $self->error_message("Failed to execute fast tier command(err: $snv_err) with params:\n" . Data::Dumper::Dumper(\%params));
        }
    }else{
        $self->debug_message("No detected variants for $name, skipping tiering");
        map {Genome::Sys->copy_file($path_to_tier, $_)}($tier1_path, $tier2_path, $tier3_path, $tier4_path);
    }
    unless(-e "$tier1_path" and -e "$tier2_path" and -e "$tier3_path" and -e "$tier4_path"){
        die $self->error_message("SNV fast tier output not found with params:\n" . (%params?(Data::Dumper::Dumper(\%params)):''));
    }
}


sub params_for_result {
    my $self = shift;
    my $qual = shift;
    my $variant_type = shift;
    my $build = $self->build;

    my @results = $build->results;
    my $target_class;
    if($qual eq 'hq') {
        $target_class = 'Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered';
    } elsif($qual eq 'lq') {
        $target_class = 'Genome::Model::Tools::DetectVariants2::Result::Combine::LqUnion';
    } else {
        die('Unknown quality passed to params_for_result: ' . $qual);
    }
    my @dv2_results = grep($_->isa($target_class), @results);
    my ($result, @extra) = grep($_->variant_type eq $variant_type, @dv2_results);
    return unless $result;
    if(@extra) {
        die $self->error_message('Multiple results unexpected for ' . $qual . ' ' . $variant_type);
    }

    return (
        variant_type => $variant_type,
        prior_result_id => $result->id,
        annotation_build_id => $build->annotation_build->id,
        classifier_version => $build->tiering_version,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    $result->add_user(label => 'uses', user => $build);

    my $effects_dir = $build->data_directory."/effects";
    unless(-d $effects_dir){
        Genome::Sys->create_directory($effects_dir);
    }

    for my $f (glob($result->output_dir . '/*')) {
        my $name = File::Basename::fileparse($f);

        Genome::Sys->create_symlink($f, join('/', $effects_dir, $name));
    }

    return 1;
}

1;

