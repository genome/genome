package Genome::Model::Tools::DetectVariants2::Result::Combine::UnionCnv;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Combine::UnionCnv{
    is => 'Genome::Model::Tools::DetectVariants2::Result::Combine',
    doc => 'Union snvs into one file',
};

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'union-cnv' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };
sub allocation_subdir_prefix { 'union_cnv' };
sub _variant_type { 'cnvs' };

sub _combine_variants {
    my $self = shift;
    # TODO Figure out how to eliminate this combine step in cnv.
    #   For now, this simply copies the cnvs.hq file from input_a into the output dir.

    my $input_a = $self->input_directory_a."/cnvs.hq";
    my $input_b = $self->input_directory_b."/cnvs.hq";

    my $a = -e $input_a;
    my $b = $self->line_count( $input_b );
    if($a and not $b){
        Genome::Sys->copy_file($input_a, $self->temp_staging_directory."/cnvs.hq");
    }
    else {
        die $self->error_message("Cnv Union operation found two cnv files, but this module currently only passes one forward.");
    }
    $self->status_message("Completed copying cnvs.hq file into output directory.");
    return 1;
}

sub _generate_standard_files {
    return 1;
}

sub _validate_output {
    my $self = shift;
    return 1;
}

1;
