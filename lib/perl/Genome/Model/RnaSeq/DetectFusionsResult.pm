package Genome::Model::RnaSeq::DetectFusionsResult;

use strict;
use warnings;

use Genome;
use File::Which;



class Genome::Model::RnaSeq::DetectFusionsResult {
    is => "Genome::SoftwareResult::Stageable",
    is_abstract => 1,
    has_param => [
        version => {
            is => 'Text',
            doc => 'the version of the fusion detector to run'
        },
        detector_params => {
            is => 'Text',
            doc => 'parameters for the chosen fusion detector'
        },
        picard_version => {
            is => 'Text',
            doc => 'the version of picard used to manipulate BAM files',
        },
    ],
    has_input => [
        alignment_result => {
            is => "Genome::SoftwareResult",
            doc => "software result from the alignment",
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'annotation build (for gene models)',
        },
    ],
};

sub _prepare_staging_directory {
    my $self = shift;

    return $self->temp_staging_directory if ($self->temp_staging_directory);

    my $tempdir = Genome::Sys->create_directory($self->output_dir . "/tmp");
    unless($tempdir) {
        die "failed to create a temp staging directory for completed files";
    }
    $self->temp_staging_directory($tempdir);

    return 1;
}

sub resolve_allocation_subdirectory {
    die "Must define resolve_allocation_subdirectory in your subclass of Genome::SoftwareResult::Stageable";
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub _remove_staging_directory {
    my $self = shift;

    Genome::Sys->remove_directory_tree($self->temp_staging_directory);

    if (-d $self->temp_staging_directory) {
        die ("You tried to remove the temp directory but it still exists...");
    }

    return 1;
}


sub _put_bowtie_version_in_path {
    my $self = shift;

    my $bowtie_path = Genome::Model::Tools::Bowtie->path_for_bowtie_version($self->alignment_result->bowtie_version);
    unless($bowtie_path){
        die($self->error_message("unable to find a path for the bowtie version specified!"));
    }

    $bowtie_path =~ s/bowtie$//;
    $ENV{PATH} = $bowtie_path . ":" . $ENV{PATH};
}

sub _path_for_command {
    my ($self, $version, $command) = @_;
    my $path = File::Which::which($command . $version);
    die ($self->error_message("You requested version ($version) of $command, but it is not supported")) unless $path;
    return $path;
}

1;
