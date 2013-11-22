package Genome::Model::Tools::DetectVariants2::Result::Combine::UnionIndel;


use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Combine::UnionIndel{
    is => 'Genome::Model::Tools::DetectVariants2::Result::Combine',
};

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'union-indel' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };
sub allocation_subdir_prefix { 'union_indel' };
sub _variant_type { 'indels' };

sub _combine_variants {
    my $self = shift;
    
    my $indels_a = $self->input_directory_a."/indels.hq.bed";
    my $indels_b = $self->input_directory_b."/indels.hq.bed";
    my $output_file = $self->temp_staging_directory."/indels.hq.bed";

    my @input_files = ($indels_a, $indels_b);

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
    my $lq_file = $self->temp_staging_directory."/indels.lq.bed";
    `touch $lq_file`;
    return 1;
}

1;
