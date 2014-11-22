package Genome::Model::ReferenceAlignment::Command::VcfSymlinks;

use strict;
use warnings;

use Genome;
use File::Basename;
use Genome::Utility::Text qw( sanitize_string_for_filesystem );

class Genome::Model::ReferenceAlignment::Command::VcfSymlinks {
    is => 'Command::V2',
    doc => "Create a directory of symlinks to builds' VCF files, suitably renamed for data transfer",
    has => [
        output_directory => {
            is=> 'String',
            doc => 'a directory to place the symlnks to SNV and Indel VCFs.',
        }, 
        models => {
            is => 'Genome::Model::ReferenceAlignment',
            is_many => 1,
            shell_args_position => 1,
            doc => 'List of models to use if searching on a list of models.'
        },
    ],
    has_optional => [
        exclude => {
            is=>'String',
            doc=>'Don\'t include models that contain this string in the name.',
            example_values => ['Pooled_Library'],
        }
    ],
};


sub help_detail {
    return "Create a directory of symlinks to builds' SNV and Indel VCF files, suitably renamed for transfer. Note that subject names will be used as the file names and may be sanitized for the filesystem.";
}


sub execute {
    my $self = shift;

    my @builds = map { $_->last_succeeded_build } $self->models;

    my $exclude = $self->exclude;
    @builds = grep { $_->model->name !~ /\Q$exclude\E/} @builds if $exclude;


    my $output_dir = $self->output_directory;
    unless(-d $output_dir) {
        Carp::croak("$output_dir is not a directory");
    }

    #now we should have all the builds
    #for each build
    #symlink each indel and SNV file into the directory.

    for my $build (@builds) {
        my $subject_name = $build->model->subject_name;
        
        #This should clean unsavory characters out of the sample name. 
        #In theory could alter the sample name to match another sample in the same cohort.
        #If this happens we will get a failed symlink.
        my $sanitized_name = sanitize_string_for_filesystem($subject_name); 
        unless($sanitized_name eq $subject_name) {
            $self->status_message("Sanitized subject name $subject_name to $sanitized_name for use in filename.");
        }
        
        if(defined $build->model->indel_detection_strategy) {
            $self->_create_symlinked_vcf_file($build, $sanitized_name, 'indels', $output_dir);
        }
        
        if(defined $build->model->snv_detection_strategy) {
            $self->_create_symlinked_vcf_file($build, $sanitized_name, 'snvs', $output_dir);
        }

    }

    return 1;
}

sub _create_symlinked_vcf_file {
    my ($self, $build, $name, $type, $output_dir) = @_;
    unless($type eq 'snvs' || $type eq 'indels') {
        die "Invalid variant type $type requested for symlink.";
    }
    my $vcf = $type eq 'snvs' ? $build->get_snvs_vcf : $build->get_indels_vcf;
    eval {
        Genome::Sys->create_symlink($vcf, "$output_dir/$name.$type.vcf.gz");
    };
    if($@) {
        die "Unable to create a symlink for the $type VCF in build " . $build->id . " for " . $build->model->subject_name . ". Check that the subject name is unique and the $type.vcf.gz file exists in the build.\n"; 
    }
    return 1;
}




1;

