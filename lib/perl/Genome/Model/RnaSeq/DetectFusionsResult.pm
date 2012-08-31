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
            is_optional => 1,
            doc => 'the version of the fusion detector to run'
        },
        detector_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'parameters for the chosen fusion detector'
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

sub resolve_allocation_subdirectory {
    die "Must define resolve_allocation_subdirectory in your subclass of Genome::SoftwareResult::Stageable";
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

#sub create {
#    my $self = shift;
#    if ($class eq __PACKAGE__ || $class->__meta__->is_abstract) {
#        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
#        return $class->SUPER::create(@_);
#    }
#
#    $self->SUPER::create(@_);
#
#    return $self;
#}

sub _put_bowtie_version_in_path {
    my $self = shift;

    my $bowtie_path = Genome::Model::Tools::Bowtie->path_for_bowtie_version($self->alignment_result->bowtie_version);
    unless($bowtie_path){
        die($self->error_message("unable to find a path for the bowtie version specified!"));
    }

    $bowtie_path =~ s/bowtie$//;
    $ENV{PATH} = $bowtie_path . ":" . $ENV{PATH};
}

sub _get_fastq_files_for_model {
    my $self = shift;

    my $alignment_result = $self->alignment_result;
    unless($alignment_result){
        die("Could not find an alignment result for build: " . $self->build->__display_name__);
    }

    $alignment_result->add_user(label => 'uses' , user => $self);

    my $fastq1 = $self->temp_staging_directory . "/fastq1";
    my $fastq2 = $self->temp_staging_directory . "/fastq2";

    my $conversion_cmd_first_read = Genome::Model::Tools::Sam::BamToFastq->create(
        bam_file  => $alignment_result->bam_file,
        fastq_file  => $fastq1,
        include_flag => 0x0040,
    );

    unless($conversion_cmd_first_read->execute){
        die("Could not convert the bam file to fastq! (first reads)");
    }

    my $conversion_cmd_second_read = Genome::Model::Tools::Sam::BamToFastq->create(
        bam_file  => $alignment_result->bam_file,
        fastq_file  => $fastq2,
        include_flag => 0x0080,
    );

    unless($conversion_cmd_second_read->execute){
        die("Could not convert the bam file to fastq! (second reads)");
    }

    return ($fastq1, $fastq2);
}

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

sub _remove_staging_directory {
    my $self = shift;

    Genome::Sys->remove_directory_tree($self->temp_staging_directory);

    if (-d $self->temp_staging_directory) {
        die ("You tried to remove the temp directory but it still exists...");
    }

    return 1;
}

sub _path_for_command {
    my ($self,$version,$command) = @_;
    my $path = File::Which::which($command . $version);
    die ($self->error_message("You requested version ($version) of $command, but it is not supported")) unless $path;
    return $path;
}

1;
