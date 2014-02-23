package Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered {
    is => 'Genome::Model::Tools::DetectVariants2::Result::Classify',
    has_input =>[
        previously_discovered_result_id => {
            is => 'Text',
            doc => 'ID of the result containing the prior variants',
        },
    ],
    has => [
        previously_discovered_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            id_by => 'previously_discovered_result_id',
        },
    ],
    has_param => [
        skip_filtering => {
            is => 'Boolean',
        },
    ],
};

sub _validate_inputs {
    my $self = shift;

    unless($self->previously_discovered_result) {
        $self->error_message('No previously discovered result found.');
        return;
    }

    unless(-e (join('/', $self->previously_discovered_result->output_dir, $self->variant_type . 's.hq.bed'))) {
        $self->error_message('Could not find ' . $self->variant_type . ' file for previously discovered variants.');
        return;
    }

    return $self->SUPER::_validate_inputs;
}

sub _classify_variants {
    my $self = shift;

    $self->previously_discovered_result->add_user(label => 'uses', user => $self);
    $self->prior_result->add_user(label => 'uses', user => $self);

    my $type = $self->variant_type;
    my $output_file = join('/', $self->temp_staging_directory, $type . 's.hq.novel.v2.bed');
    my $previously_detected_output_file = join('/', $self->temp_staging_directory, $type . 's.hq.previously_detected.v2.bed');

    my $previously_discovered_path = $self->previously_discovered_result->path($type . 's.hq.bed');
    my $prior_path = $self->prior_result->path($type . 's.hq.bed');

    unless (-s $prior_path){
    
        $self->debug_message("high confidence input is empty, skipping intersection");
        File::Copy::copy($prior_path, $output_file);
        File::Copy::copy($prior_path, $previously_detected_output_file);
        return 1;
    }
    if ($self->skip_filtering) {
        $self->debug_message("Skipping filtering");
        File::Copy::copy($prior_path, $output_file);
        my $fh = Genome::Sys->open_file_for_writing($previously_detected_output_file);
        $fh->close;
    }
    else {
        $self->debug_message("Not skipping filtering");
        my $snv_compare = Genome::Model::Tools::Joinx::Intersect->create(
            input_file_a => $prior_path,
            input_file_b => $previously_discovered_path,
            miss_a_file => $output_file,
            output_file => $previously_detected_output_file,
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
    }
    

    return 1;
}

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'dv2-previously-discovered-result' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    return join('/', 'build_merged_alignments', $self->id, 'dv2-classify-previously-discovered-' . $staged_basename);
};

sub path {
    my $self = shift;
    my ($str) = @_;

    my $type = $self->variant_type;

    if($str eq $type . 's.hq.bed') {
        return $self->path($type . 's.hq.novel.v2.bed');
    } elsif($str eq $type . 's.lq.bed') {
        return $self->path($type . 's.hq.previously_detected.v2.bed');
    } else {
        return $self->SUPER::path(@_);
    }
}
