package Genome::Model::Tools::DetectVariants2::Result::Combine::UnionSnv;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Combine::UnionSnv{
    is => 'Genome::Model::Tools::DetectVariants2::Result::Combine',
    doc => 'Union snvs into one file',
};

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'union-snv' };
sub resolve_allocation_disk_group_name { 'info_genome_models' };
sub allocation_subdir_prefix { 'union_snv' };
sub _variant_type { 'snvs' };

sub _combine_variants {
    my $self = shift;
    my $snvs_a = $self->input_directory_a."/snvs.hq.bed";
    my $snvs_b = $self->input_directory_b."/snvs.hq.bed";
    my $output_file = $self->temp_staging_directory."/snvs.hq.bed";

    my @input_files = ($snvs_a, $snvs_b);

    # Using joinx with --merge-only will do a union, effectively
    my $union_command = Genome::Model::Tools::Joinx::Sort->create(
        input_files => \@input_files,
        merge_only => 1,
        output_file => $output_file,
    );
    
    unless ($union_command->execute) {
        $self->error_message("Error executing union command");
        die $self->error_message;
    }

    # When unioning, there is no "fail" really, everything should be in the hq file
    my $lq_file = $self->temp_staging_directory."/snvs.lq.bed";
    `touch $lq_file`;
    
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

    my $merge_cmd = Genome::Model::Tools::Joinx::VcfMerge->create(
        input_files => [ ($input_a_vcf,$input_b_vcf)],
        output_file => $output_file,
        merge_samples => 1,
        clear_filters => 1,
        use_bgzip => 1,
        use_version => "1.6",
    );
    unless($merge_cmd->execute){
        die $self->error_message("Could not complete merge operation.");
    }

    return 1;
}

1;
