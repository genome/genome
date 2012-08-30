package Genome::Model::Tools::DetectVariants2::Result::Combine::IntersectSnv;

use warnings;
use strict;

use Genome;


class Genome::Model::Tools::DetectVariants2::Result::Combine::IntersectSnv{
    is => 'Genome::Model::Tools::DetectVariants2::Result::Combine',
    doc => 'intereset snv filter results into one file',
};

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'intersect-snv' };
sub resolve_allocation_disk_group_name { 'info_genome_models' };
sub allocation_subdir_prefix { 'intersect_snv' };
sub _variant_type { 'snvs' };

sub _combine_variants {
    my $self = shift;

    my $snvs_a = $self->input_directory_a."/snvs.hq.bed";
    my $snvs_b = $self->input_directory_b."/snvs.hq.bed";
    my $output_file = $self->temp_staging_directory."/snvs.hq.bed";
    my $miss_a_file = $self->temp_staging_directory."/snvs.lq.a.bed";
    my $miss_b_file = $self->temp_staging_directory."/snvs.lq.b.bed";

    my $scratch_dir = Genome::Sys->create_temp_directory;

    my $sorted_snvs_a = $scratch_dir."/a.hq.bed";
    my $sorted_snvs_b = $scratch_dir."/b.hq.bed";

    for ([$snvs_a, $sorted_snvs_a], [$snvs_b, $sorted_snvs_b]) {
        my $sort_snvs_cmd = Genome::Model::Tools::Joinx::Sort->create(
            input_files => [$_->[0]],
            output_file => $_->[1],
        );
        unless($sort_snvs_cmd->execute){
            die $self->error_message("Error executing sort command on " . $_->[0]);
        }
    }

    my $intersect_command = Genome::Model::Tools::Joinx::Intersect->create(
        input_file_a => $sorted_snvs_a,
        input_file_b => $sorted_snvs_b,
        output_file => $output_file,
        miss_a_file => $miss_a_file,
        miss_b_file => $miss_b_file,
    );

    unless ($intersect_command->execute) {
        $self->error_message("Error executing intersect command");
        die $self->error_message;
    }

    # Create an "lq" file that has things that were either in only file a or only file b
    # Using joinx with --merge-only will do a union, effectively
    my $lq_file = $self->temp_staging_directory."/snvs.lq.bed";
    my $merge_cmd = Genome::Model::Tools::Joinx::Sort->create(
        merge_only => 1,
        input_files => [$miss_a_file, $miss_b_file],
        output_file => $lq_file,
    );
    unless ($merge_cmd->execute) {
        $self->error_message("Failed to combine $miss_a_file and $miss_b_file into $lq_file");
        die $self->error_message;
    }

    $self->_generate_vcf;

    return 1;
}

sub _generate_vcf {
    my $self = shift;
    my $input_a_vcf = $self->input_directory_a."/snvs.vcf.gz";
    unless(-s $input_a_vcf){
        $self->status_message("Could not find vcf at: ".$input_a_vcf." not creating a vcf for this operation.");
        return;
    }
    my $input_b_vcf = $self->input_directory_b."/snvs.vcf.gz";
    unless(-s $input_b_vcf){
        $self->status_message("Could not find vcf at: ".$input_b_vcf." not creating a vcf for this operation.");
        return;
    }
    my $output_file = $self->temp_staging_directory."/snvs.vcf.gz";
    my $vcf_files;
    my $vcf_a_source = $self->get_vcf_source($input_a_vcf);
    my $vcf_b_source = $self->get_vcf_source($input_b_vcf);
    my $source_ids;

    # FIXME In the future, probably just take the order in the processing profile
    # FIXME but for now... this was done in germline. But maybe in somatic we want to prioritize sniper over samtools? 
    if($vcf_a_source =~ m/samtools/i){
        #if A has samtools, then the ordering is fine, change nothing
    } elsif ( $vcf_b_source =~ m/samtools/i) {
        #if B has samtools, swap the ordering
        ($input_a_vcf,$input_b_vcf) = ($input_b_vcf,$input_a_vcf);
        ($vcf_a_source,$vcf_b_source) = ($vcf_b_source,$vcf_a_source);
    } else {
        #if we cannot locate samtools, die
        die $self->error_message("Could not positively identify samtools input!");
    }
    my $merge_cmd = Genome::Model::Tools::Joinx::VcfMerge->create(
        input_files => [ ($input_a_vcf,$input_b_vcf)],
        output_file => $output_file,
        merge_samples => 1,
        clear_filters => 1,
        use_bgzip => 1,
        use_version => "1.6",
        ratio_filter => "1.0,IntersectionFailure,Variant callers do not agree on this position",
    );

    unless($merge_cmd->execute){
        die $self->error_message("Could not complete merge operation.");
    }

    return 1;
}

sub get_vcf_source {
    my $self = shift;
    my $vcf = shift;

    unless(-s $vcf){
        die $self->error_message("Cannot determine the source of a file that doesn't exist...");
    }

    my $fh = Genome::Sys->open_gzip_file_for_reading($vcf);
    my $source;
    while(my $line = $fh->getline){
        chomp $line;
        if( ($line =~ m/\#\#/) && ($line =~ m/source/)) {
            (undef,$source) = split /\=/, $line;
            last;
        }
        unless( ($line =~ m/\#\#/) ){
            die $self->error_message("Could not find source tag in the header.");
        }
    }
    return $source;
}

sub _validate_output {
    my $self = shift;
    my $variant_type = $self->_variant_type;
    my $input_a_file = $self->input_directory_a."/".$variant_type.".hq.bed";
    my $input_b_file = $self->input_directory_b."/".$variant_type.".hq.bed";
    my $hq_output_file = $self->temp_staging_directory."/".$variant_type.".hq.bed";
    my $lq_output_file = $self->temp_staging_directory."/".$variant_type.".lq.bed";
    my $input_total = $self->line_count($input_a_file) + $self->line_count($input_b_file);

    # Count hq * 2 because every hq line for an intersection implies 2 lines from input combined into one
    my $output_total = $self->line_count($hq_output_file) * 2 + $self->line_count($lq_output_file);
    unless(($input_total - $output_total) == 0){
        die $self->error_message("Combine operation in/out check failed. Input total: $input_total \toutput total: $output_total");
    }
    return 1;
}

1;
