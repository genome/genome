package Genome::Model::SomaticVariation::Command::IdentifyPreviouslyDiscoveredVariations;

use strict;
use warnings;
use Genome;
use Genome::Model::Tools::DetectVariants2::Utilities qw(
    final_result_for_variant_type
);

class Genome::Model::SomaticVariation::Command::IdentifyPreviouslyDiscoveredVariations{
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
        }
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
    ],
};

sub shortcut {
    my $self = shift;

    if($self->should_skip_run) {
        for my $type ('snv', 'indel') {
            $self->skip_run($type);
        }
        return 1;
    }

    my $build = $self->build;

    my @results;
    my $expected_result_count = 0;

    for my $type ('snv', 'indel') {
        my $accessor = join('_', $type, 'detection_strategy');
        if($build->model->$accessor) {
            $expected_result_count++;
            my @params_for_result = $self->params_for_result($type);
            return unless @params_for_result; #can't use a result (typically previous data was not in a result)
            my $result = Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered->get_with_lock(@params_for_result);
            push @results, $result if $result;
        }
    }

    unless($expected_result_count > 0 && scalar(@results) == $expected_result_count) {
        return;
    }

    for my $r (@results) {
        $self->link_result_to_build(@results) or die('Failed to link result.');
    }

    return 1;
}

sub execute {
    my $self  = shift;

    if($self->should_skip_run) {
        for my $type ('snv', 'indel') {
            $self->skip_run($type);
        }
        return 1;
    }

    my $build = $self->build;

    $self->status_message("Comparing detected variants to previously discovered variations");

    my ($snv_result, $indel_result);

    my $prev_variations_build = $build->previously_discovered_variations_build;
    $snv_result   = $prev_variations_build->snv_result;
    $indel_result = $prev_variations_build->indel_result;

    unless ($indel_result or $snv_result) {
        die $self->error_message("No indel or snv result found on previously discovered variations build. This is unsupported!  Failing.");
    }

    my $version = 2;
    #my $version = GMT:BED:CONVERT::version();  TODO, something like this instead of hardcoding

    if ($build->snv_detection_strategy){
        my @params = $self->params_for_result('snv');
        $self->status_message(Data::Dumper::Dumper(\@params));
        if(@params) {
            my $result = Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered->get_or_create(@params);
            $self->status_message("Using result from get_or_create ".$result->id);
            $self->link_result_to_build($result);
        } elsif($snv_result) {

            my $detected_snv_path = defined($build->loh_version) ? $build->data_set_path("loh/snvs.somatic",$version,'bed') : $build->data_set_path("variants/snvs.hq",$version,"bed"); 
            my $novel_detected_snv_path = $build->data_set_path("novel/snvs.hq.novel",$version,'bed');
            my $previously_detected_snv_path = $build->data_set_path("novel/snvs.hq.previously_detected",$version,'bed');

            my $snv_result_path = join('/', $snv_result->output_dir, 'snvs.hq.bed');

            unless (-e $snv_result_path){
                die $self->error_message("Snv feature list does not have an associated file!");
            }

            unless (-e $detected_snv_path){
                die $self->error_message("No high confidence detected snvs to filter against previously discovered variants");
            }

            if ($build->processing_profile->filter_previously_discovered_variants) {
                if (-s $detected_snv_path){
                    my $snv_output_tmp_file = Genome::Sys->create_temp_file_path();
                    my $previously_detected_output_tmp_file = Genome::Sys->create_temp_file_path();
                    my $snv_compare = Genome::Model::Tools::Joinx::Intersect->create(
                        input_file_a => $detected_snv_path,
                        input_file_b => $snv_result_path,
                        miss_a_file => $snv_output_tmp_file,
                        output_file => $previously_detected_output_tmp_file,
                        dbsnp_match => 1,
                    );
                    unless ($snv_compare){
                        die $self->error_message("Couldn't create snv comparison tool!");
                    }
                    my $snv_rv = $snv_compare->execute();
                    my $snv_err = $@;
                    unless ($snv_rv){
                        die $self->error_message("Failed to execute snv comparison(err: $snv_err )");
                    }
                    $self->status_message("Intersection against previously discovered snv feature list complete");
                    File::Copy::copy($snv_output_tmp_file, $novel_detected_snv_path);
                    File::Copy::copy($previously_detected_output_tmp_file, $previously_detected_snv_path);
                }
                else{
                    $self->status_message("high confidence snv output is empty, skipping intersection");
                    Genome::Sys->create_directory($build->data_directory."/novel");
                    File::Copy::copy($detected_snv_path, $novel_detected_snv_path);
                    File::Copy::copy($detected_snv_path, $previously_detected_snv_path);
                }
            }
            else {
                $self->status_message("Skipping filtering");
                Genome::Sys->create_directory($build->data_directory."/novel");
                File::Copy::copy($detected_snv_path, $novel_detected_snv_path);
                my $fh = Genome::Sys->open_file_for_writing($previously_detected_snv_path);
                $fh->close;
            }
        }
        else{
            $self->status_message("No snv feature list found on previously discovered variations build, skipping snv intersection");
            $self->skip_run('snv');
        }
    }

    if ($build->indel_detection_strategy) {
        if(my @params = $self->params_for_result('indel')) {
            my $result = Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered->get_or_create(@params);
            $self->link_result_to_build($result);
        } elsif($indel_result) {

            my $detected_indel_path            = $build->data_set_path("variants/indels.hq",$version,"bed"); 
            my $novel_detected_indel_path      = $build->data_set_path("novel/indels.hq.novel",$version,"bed");
            my $previously_detected_indel_path = $build->data_set_path("novel/indels.hq.previously_detected", $version, "bed");

            my $indel_result_path = join('/', $indel_result->output_dir, 'indels.hq.bed');

            unless (-e $indel_result_path){
                die $self->error_message("Indel feature list does not have an associated file!");
            }

            unless (-e $detected_indel_path){
                die $self->error_message("No high confidence detected indels to filter against previously discovered variants");
            }

            if (-s $detected_indel_path){
                my $indel_output_tmp_file = Genome::Sys->create_temp_file_path();
                my $previously_detected_output_tmp_file = Genome::Sys->create_temp_file_path();
                my $indel_compare = Genome::Model::Tools::Joinx::Intersect->create(
                    input_file_a => $detected_indel_path,
                    input_file_b => $indel_result_path,
                    miss_a_file => $indel_output_tmp_file,
                    output_file => $previously_detected_output_tmp_file,
                    exact_allele => 1,
                );
                unless ($indel_compare){
                    die $self->error_message("Couldn't create indel comparison tool!");
                }
                my $indel_rv = $indel_compare->execute();
                my $indel_err = $@;
                unless ($indel_rv){
                    die $self->error_message("failed to execute indel comparison(err: $indel_err )");
                }
                $self->status_message("intersection against previously discovered indel feature list complete");
                File::Copy::copy($indel_output_tmp_file, $novel_detected_indel_path);
                File::Copy::copy($previously_detected_output_tmp_file, $previously_detected_indel_path);
            }
            else{
                $self->status_message("high confidence indel output is empty, skipping intersection");
                Genome::Sys->create_directory($build->data_directory."/novel");
                File::Copy::copy($detected_indel_path, $novel_detected_indel_path);
                File::Copy::copy($detected_indel_path, $previously_detected_indel_path);
            }
        }
        else{
            $self->status_message("No indel feature list found on previously discovered variations build, skipping indel intersection");
            $self->skip_run('indel');
        }
    }

    $self->status_message("Identify Previously Discovered Variations step completed");
    return 1;
}

sub should_skip_run {
    my $self = shift;
    my $build = $self->build;

    unless ($build){
        die $self->error_message("no build provided!");
    }

    unless(defined($build->model->snv_detection_strategy) or defined($build->model->indel_detection_strategy)){
        $self->status_message("No SNV or indel detection strategy, skipping identify previously discovered variants.");
        return 1;
    }

    my $prev_variations_build = $build->previously_discovered_variations_build;
    unless ($prev_variations_build) {
        $self->warning_message('no previously_discovered_variations_build provided !');
        return 1;
    }

    return;
}

sub skip_run {
    my $self = shift;
    my $variant_type = shift;
    my $build = $self->build;

    Genome::Sys->create_directory($build->data_directory."/novel");

    my $accessor = join('_', $variant_type, 'detection_strategy');
    return 1 unless $build->model->$accessor; #don't make empty files if we're not running this

    my $novel = join('/', $build->data_directory, 'novel', $variant_type . 's.hq.novel.v2.bed');
    my $previously_detected = join('/', $build->data_directory, 'novel', $variant_type . 's.hq.previously_detected.v2.bed');

    my $detected_path;
    if(my %params = $self->params_for_result($variant_type)) {
        my $r = Genome::Model::Tools::DetectVariants2::Result::Base->get($params{prior_result_id});
        $detected_path = $r->path($variant_type . '.hq.bed');
    } else {
        my $version = '2';
        $detected_path = ($variant_type eq 'snv' && defined($build->loh_version)) ? $build->data_set_path("loh/snvs.somatic",$version,'bed') : $build->data_set_path("variants/" . $variant_type . "s.hq",$version,"bed");
    }
    File::Copy::copy($detected_path, $novel);
    system("touch $previously_detected");

    return 1;
}


sub params_for_result {
    my $self = shift;
    my $variant_type = shift;
    my $build = $self->build;

    my $prior_result;
    if($variant_type eq 'snv' and $build->loh_version) {
        my @results = $build->results;
        for my $r (@results) {
            if($r->class eq 'Genome::Model::Tools::DetectVariants2::Classify::Loh') {
                $prior_result = $r;
                last;
            }
        }
    } else {
        $prior_result = final_result_for_variant_type([$build->results], $variant_type);
    }

    my $variations_build = $build->previously_discovered_variations_build;
    return unless $variations_build;

    my $accessor = $variant_type . '_result';
    my $previously_discovered_result = $variations_build->$accessor;

    unless($prior_result and $previously_discovered_result) {
        return; #can't create a result for old things
    }

    my $skip = 1;
    if ($build->processing_profile->filter_previously_discovered_variants) {
        $skip = 0;
    }

    return (
        prior_result_id => $prior_result->id,
        previously_discovered_result_id => $previously_discovered_result->id,
        classifier_version => 1,
        variant_type => $variant_type,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        skip_filtering => $skip,
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    $self->status_message('Linking result ' . $result->id . ' to build.');
    $result->add_user(user => $build, label => 'uses');
    Genome::Sys->create_directory($build->data_directory."/novel");

    my $type = $result->variant_type;
    my $novel = $type . 's.hq.novel.v2.bed';
    my $previously_detected = $type . 's.hq.previously_detected.v2.bed';
    my $novel_detected_path      = join('/', $build->data_directory, 'novel', $novel);
    my $previously_detected_path = join('/', $build->data_directory, 'novel', $previously_detected);
    Genome::Sys->create_symlink($result->path($novel), $novel_detected_path);
    Genome::Sys->create_symlink($result->path($previously_detected), $previously_detected_path);

    return 1;
}

1;

