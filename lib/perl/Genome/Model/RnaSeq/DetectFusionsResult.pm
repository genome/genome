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
        picard_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'the version of picard used to manipulate BAM files',
            default_value => '1.82',
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
    has => [
        _fastq_files => {
            # returns an array ref of [fastq1, fastq2]
            is_constant => 1,
            calculate => q($self->_get_fastq_files();),
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

sub _get_fastq_files {
    my $self = shift;

    my $alignment_result = $self->alignment_result;
    unless($alignment_result){
        die("Could not find an alignment result for build: " . $self->build->__display_name__);
    }

    $alignment_result->add_user(label => 'uses' , user => $self);

    
    my $fastq1 = $self->temp_staging_directory . "/fastq1";
    my $fastq2 = $self->temp_staging_directory . "/fastq2";
    
    my $cmd = Genome::Model::Tools::Picard::StandardSamToFastq->create(
        input            => $alignment_result->bam_file,
        fastq            => $fastq1,
        second_end_fastq => $fastq2,
        re_reverse       => 1,
        include_non_pf_reads => 1,
        include_non_primary_alignments => 0,
        use_version      => $self->picard_version,
        maximum_memory => 30,
        maximum_permgen_memory => 256,
        max_records_in_ram => 5_000_000,
        additional_jvm_options => '-XX:-UseGCOverheadLimit',
    );
    unless ($cmd->execute){
        die("Could not convert the bam file to fastq!");
    }
    my @result = ($fastq1, $fastq2);
    return \@result;
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
