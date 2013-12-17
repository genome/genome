package Genome::Model::Event::Build::ReferenceAlignment::LaneQc::CompareSnps;

use strict;
use warnings;

use Genome;
require File::Path;
use Cwd;

class Genome::Model::Event::Build::ReferenceAlignment::LaneQc::CompareSnps {
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
        # TODO why isn't this a die or a return?
        $self->status_message("No valid gold_snp_path for the build, aborting compare SNPs!");
    }

    my $output_dir = $build->qc_directory;
    File::Path::mkpath($output_dir) unless (-d $output_dir);
    unless (-d $output_dir) {
        die $self->error_message("Failed to create output_dir ($output_dir).");
    }

    my $geno_path = $self->resolve_geno_path_for_build($build);

    unless (-s $geno_path) {
        $self->warning_message("Genotype file is empty: $geno_path. Joinx intersect did not find any intersection.");
        return 1;
    }

    #TODO: Remove Over-Ambiguous Glob
    my @variant_files = glob($build->variants_directory . '/snv/samtools-*/snvs.hq');
    unless(scalar @variant_files eq 1) {
        die $self->error_message("Could not find samtools output for run.");
    }
    my $variant_file = $variant_files[0];
    unless ( -s $variant_file ) {
        die $self->error_message("Variant file missing/empty: $variant_file");
    }
    $variant_file = Cwd::abs_path($variant_file);

    my %compare_snps_result_params = (
        genotype_file => $geno_path,
        variant_file => $variant_file,
        sample_name => $model->subject->name,
        reference_build => $build->reference_sequence_build->version,
    );
    if ($build->region_of_interest_set_name) {
        $compare_snps_result_params{bam_file} = $build->whole_rmdup_bam_file;
        $compare_snps_result_params{flip_alleles} = 1;
    }
    my $result = Genome::Model::Tools::Analysis::LaneQc::CompareSnpsResult->get_or_create(%compare_snps_result_params);
    unless ($result) {
        die $self->error_message("Failed to create Genome::Model::Tools::Analysis::LaneQc::CompareSnpsResult command.");
    }

    $result->add_user( user_id => $build->id, user_class_name => $build->class, label => 'uses' );

    die 'Missing args for creating symlink' unless $result->output_file and $build->compare_snps_file;
    Genome::Sys->create_symlink_and_log_change($self, $result->output_file, $build->compare_snps_file);

    my $metrics_rv = Genome::Model::ReferenceAlignment::Command::CreateMetrics::CompareSnps->execute(
        build_id => $self->build_id,
    );
    Carp::confess "Could not create compare_snps metrics for build " . $self->build_id unless $metrics_rv;

    return 1;
}

sub resolve_geno_path_for_build {
    my $self = shift;
    my $build = shift;

    my $geno_path;
    if ($build->region_of_interest_set_name) {
        my $output_dir = $build->qc_directory;

        my $feature_list = Genome::FeatureList->get(name => $build->region_of_interest_set_name);
        # This cleans up the FeatureList, e.g. removes chr prefix.
        my $roi_bed_path = $feature_list->processed_bed_file(short_name => 0);

        # Convert original gold2geno file into a BED for easy intersection with FeatureList
        my $gold2geno_bed_path = UR::Value::FilePath->get("$output_dir/genotype.gold2geno.bed");
        my $original_gold2geno_path = UR::Value::FilePath->get($build->gold_snp_build->gold2geno_file_path);
        my $gold2geno_bed = Genome::Sys->open_file_for_writing($gold2geno_bed_path);
        my $gold2geno = Genome::Sys->open_file_for_reading($original_gold2geno_path);
        while (my $line = $gold2geno->getline) {
            chomp($line);
            my ($chrom, $end, $allele) = split("\t", $line);
            my $start = $end - 1;
            $gold2geno_bed->print(join("\t", $chrom, $start, $end, $allele), "\n");
        }
        $gold2geno_bed->close;
        $gold2geno->close;
        unless ($original_gold2geno_path->line_count == $gold2geno_bed_path->line_count) {
            die "Converted gold2geno BED ($gold2geno_bed_path) line count does not match original gold2geno line count ($original_gold2geno_path)";
        }

        my $sorted_bed_path = Genome::Sys->create_temp_file_path();

        my $sort_cmd = Genome::Model::Tools::Joinx::Sort->create(
            input_files => [$roi_bed_path],
            output_file => "$sorted_bed_path",
        );

        unless ($sort_cmd->execute) {
            die $self->error_message("Failed to sort feature list BED.");
        }

        my $sorted_gold2geno_bed_path = Genome::Sys->create_temp_file_path;
        my $rv = Genome::Model::Tools::Joinx::Sort->execute(
            input_files => [$gold2geno_bed_path],
            output_file => "$sorted_gold2geno_bed_path",
        );
        unless ($rv) {
            $self->error_message("Failed to sort gold2geno bed");
            return;
        }

        my $intersected_gold2geno_path = UR::Value::FilePath->get("$output_dir/intersected_genotype.gold2geno");
        my $intersected_gold2geno_bed_path = UR::Value::FilePath->get("$intersected_gold2geno_path.bed");
        my $intersect_cmd = Genome::Model::Tools::Joinx::Intersect->create(
            input_file_a => "$sorted_gold2geno_bed_path", # genotype first
            input_file_b => "$sorted_bed_path",
            output_file  => "$intersected_gold2geno_bed_path",
        );
        unless ($intersect_cmd->execute) {
            die $self->error_message("Failed to intersect sorted feature list and genotype BEDs.");
        }

        # Convert the intersected gold2geno BED back to a "standard" gold2geno file
        my $intersected_gold2geno_bed = Genome::Sys->open_file_for_reading($intersected_gold2geno_bed_path);
        my $intersected_gold2geno = Genome::Sys->open_file_for_writing($intersected_gold2geno_path);
        while (my $line = $intersected_gold2geno_bed->getline) {
            chomp($line);
            my ($chrom, $start, $end, $allele) = split("\t", $line);
            $intersected_gold2geno->print(join("\t", $chrom, $end, $allele), "\n");
        }
        $intersected_gold2geno_bed->close;
        $intersected_gold2geno->close;
        unless ($intersected_gold2geno_path->line_count == $intersected_gold2geno_bed_path->line_count) {
            die "Converted gold2geno ($intersected_gold2geno_path) line count does not match BED line count ($intersected_gold2geno_bed_path)";
        }
        $geno_path = "$intersected_gold2geno_path";
    } else {
        $geno_path = $build->gold_snp_build->gold2geno_file_path;
    }

    unless ( -e $geno_path ) {
        die $self->error_message("Genotype file missing: $geno_path");
    }

    return $geno_path;
}

sub validate_gold_snp_path {
    my $self = shift;

    my $gold_snp_path = $self->build->gold_snp_path;
    unless ($gold_snp_path and -s $gold_snp_path) {
        $self->status_message('No gold_snp_path provided for the build or it is empty');
        return;
    }

    my $head    = `head -1 $gold_snp_path`;
    my @columns = split /\s+/, $head;

    unless (@columns and @columns == 9) {
        $self->status_message("Gold snp file: $gold_snp_path is not 9-column format");
        return;
    }
    return 1;
}

1;
