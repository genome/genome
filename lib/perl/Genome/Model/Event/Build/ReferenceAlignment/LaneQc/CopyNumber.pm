package Genome::Model::Event::Build::ReferenceAlignment::LaneQc::CopyNumber;

use strict;
use warnings;

use Genome;
use File::Path;

class Genome::Model::Event::Build::ReferenceAlignment::LaneQc::CopyNumber {
    is => [ 'Genome::Model::Event' ],
};

sub execute {
    my $self  = shift;
    my $model = $self->model;
    my $build = $self->build;

    my @instrument_data = $build->instrument_data;
    if (@instrument_data > 1) {
        my $package = __PACKAGE__;
        die $self->error_message("Build has many instrument data, $package is designed to run on a per-lane basis.");
    }

    if ( !$self->validate_gold_snp_path ) {
        die $self->error_message("No valid gold_snp_path for the build, aborting copy number!");
    }

    my $output_dir = $build->data_directory . '/qc';
    mkpath($output_dir) unless (-d $output_dir);
    unless (-d $output_dir) {
        die $self->error_message("Failed to create output_dir ($output_dir).");
    }

    my $cmd = Genome::Model::Tools::Analysis::LaneQc::CopyNumber->create(
        build_id => $self->build->id,
        output_file_prefix => "$output_dir/copy_number_",
        lsf => 0,
    );
    unless ($cmd) {
        die $self->error_message("Failed to create Genome::Model::Tools::Analysis::LaneQc::CopyNumber command.");
    }

    my $cmd_executed = eval { $cmd->execute };
    unless ($cmd_executed) {
        if ($@) {
            die $self->error_message("Failed to execute CopyNumber QC analysis! Received error: $@");
        }
        else {
            die $self->error_message("Failed to execute CopyNumber QC analysis!");
        }
    }

    return 1;
}

sub validate_gold_snp_path {
    my $self = shift;

    my $gold_snp_path = $self->build->gold_snp_path;
    unless ($gold_snp_path) {
        $self->debug_message('No gold_snp_path provided for the build');
        return;
    }
    unless (-s $gold_snp_path) {
        $self->debug_message('gold_snp_path is empty ' . $gold_snp_path);
        return;
    }

    my $head    = `head -1 $gold_snp_path`;
    my @columns = split /\s+/, $head;

    unless (@columns and @columns == 9) {
        $self->debug_message("Gold snp file: $gold_snp_path is not 9-column format");
        return;
    }
    return 1;
}

1;
