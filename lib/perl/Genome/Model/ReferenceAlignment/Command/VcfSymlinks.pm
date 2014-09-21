package Genome::Model::ReferenceAlignment::Command::VcfSymlinks;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::ReferenceAlignment::Command::VcfSymlinks {
    is => 'Genome::Command::Base',
    doc => "Create a directory of symlinks to builds' VCF files, suitably renamed for data transfer",
    has => [
        output_directory => {
            is=> 'String',
            doc => 'a directory to place the symlnks',
        }, 
        builds => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            is_many => 1,
            shell_args_position => 1,
            doc => 'List of builds to use if searching on a list of builds'
        },
    ],
    has_optional => [
        exclude => {
            is=>'String',
            doc=>'Don\'t include models that contain this string in the name (ie. "Pooled_Library")',
        }
    ],
};


sub help_detail {
    return "Create a directory of symlinks to builds' VCF files, suitably renamed for transfer.";
}


sub execute {
    my $self = shift;

    my @builds = grep {$_ == $_->model->last_succeeded_build} $self->builds;

    my $exclude = $self->exclude;
    @builds = grep {not $_->model->name =~ /$exclude/} @builds if $exclude;


    my $output_dir = $self->output_directory;
    unless(-d $output_dir) {
        Carp::croak("$output_dir is not a directory");
    }

    #now we should have all the builds
    #for each build
    #symlink each indel and SNV file into the directory.
    #assuming that sample names are valid file names. I think we sanitize them such that this is true

    for my $build (@builds) {
        my $sample_name = $build->model->subject_name;
        my $indel_vcf = $build->data_directory . "/variants/indels.vcf.gz";
        Genome::Sys->create_symlink($indel_vcf, "$output_dir/$sample_name.indels.vcf.gz");
        my $snv_vcf = $build->data_directory . "/variants/snvs.vcf.gz";
        Genome::Sys->create_symlink($indel_vcf, "$output_dir/$sample_name.snvs.vcf.gz");
    }

    return 1;
}

1;

