package Genome::Model::Tools::DetectVariants2::Result::Combine::UnionuniqueSnv;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Combine::UnionuniqueSnv{
    is => 'Genome::Model::Tools::DetectVariants2::Result::Combine',
    doc => 'Union snvs into one file',
};

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'unionunique-snv' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };
sub allocation_subdir_prefix { 'unionunique_snv' };
sub _variant_type { 'snvs' };

sub _combine_variants {
    my $self = shift;
    my $snvs_a = $self->input_directory_a."/snvs.hq.bed";
    my $snvs_b = $self->input_directory_b."/snvs.hq.bed";
    my $output_file = $self->temp_staging_directory."/snvs.hq.bed";

    my @input_files = ($snvs_a, $snvs_b);

    # Using joinx with --merge-only will do a union, effectively
    my $union_command = Genome::Model::Tools::Joinx::Union->create(
        input_files => \@input_files,
        output_file => $output_file,
    );
    
    unless ($union_command->execute) {
        $self->error_message("Error executing union command");
        die $self->error_message;
    }

    # When unioning, there is no "fail" really, everything should be in the hq file
    my $lq_file = $self->temp_staging_directory."/snvs.lq.bed";
    `touch $lq_file`;
    return 1;
}

sub _validate_output {
    my $self = shift;
    my $variant_type = $self->_variant_type;
    my $input_a_file = $self->input_directory_a."/".$variant_type.".hq.bed";
    my $input_b_file = $self->input_directory_b."/".$variant_type.".hq.bed";
    my $hq_output_file = $self->temp_staging_directory."/".$variant_type.".hq.bed";
    my $lq_output_file = $self->temp_staging_directory."/".$variant_type.".lq.bed";

    unless (-e $input_a_file and -e $input_b_file and -e $hq_output_file and -e $lq_output_file) {
        $self->debug_message("Skipping _validate_output because not all bed files exist. This may be a vcf-only detector/filter strategy");
        return 1;
    }

    my $input_total = $self->line_count($input_a_file) + $self->line_count($input_b_file);

    # Count hq * 2 because every hq line for an intersection implies 2 lines from input combined into one
    my $output_total = $self->line_count($hq_output_file) + $self->line_count($lq_output_file);

    # Since we throw out non-unique variants in the union... we have to find out how many things were tossed out
    # We can figure that out by finding out how many things intersected.
    my $scratch_dir = File::Temp::tempdir(
        'Genome-Model-Tools-DetectVariants2-Combine-UnionuniqueSnv-XXXXX',
        DIR => Genome::Sys->base_temp_directory(),
        CLEANUP => 1
        );
    my $temp_intersect_file = $scratch_dir . "/UnionuniqueSnv.intersected";
    my $sorted_input_a_file = $scratch_dir . "/UnionuniqueSnv.a.sorted";
    my $sorted_input_b_file = $scratch_dir . "/UnionuniqueSnv.b.sorted";
    my $sort_command_a = Genome::Model::Tools::Joinx::Sort->create(
        input_files => [$input_a_file],
        output_file => $sorted_input_a_file,
    );
    unless($sort_command_a->execute()) {
        die "Failed to sort $input_a_file, $@";
    }
    my $sort_command_b = Genome::Model::Tools::Joinx::Sort->create(
        input_files => [$input_b_file],
        output_file => $sorted_input_b_file,
    );
    unless($sort_command_b->execute()) { 
        die "Failed to sort $input_b_file, $@";
    }
    my $intersect_command = Genome::Model::Tools::Joinx::Intersect->create(
        input_file_a => $sorted_input_a_file,
        input_file_b => $sorted_input_b_file,
        output_file => $temp_intersect_file,
        exact_pos => 1,
        exact_allele => 1,
    );
    unless ($intersect_command->execute) {
        die $self->error_message("Failed to execute intersect command to validate output");
    }
    my $offset_lines = $self->line_count($temp_intersect_file);

    # We also need to find out if there were any duplicates within each of the individual files
    # Using joinx with --merge-only will do a union, effectively
    my $temp_unique_a = $scratch_dir . "/file_a.unique";
    my $unique_a = Genome::Model::Tools::Joinx::Union->create(
        input_files => [$input_a_file],
        output_file => $temp_unique_a,
    );
    unless ($unique_a->execute) {
        die $self->error_message("Failed to execute unique_a command to validate output");
    }
    $offset_lines += $self->line_count($input_a_file) - $self->line_count($temp_unique_a);

    my $temp_unique_b = $scratch_dir . "/file_b.unique";
    my $unique_b = Genome::Model::Tools::Joinx::Union->create(
        input_files => [$input_b_file],
        output_file => $temp_unique_b,
    );
    unless ($unique_b->execute) {
        die $self->error_message("Failed to execute unique_b command to validate output");
    }
    $offset_lines += $self->line_count($input_b_file) - $self->line_count($temp_unique_b);

    unless(($input_total - $output_total - $offset_lines) == 0){
        die $self->error_message("Combine operation in/out check failed. Input total: $input_total \toutput total: $output_total\t with an intersected offset of $offset_lines");
    }
    return 1;
}

1;
