package Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult;

use strict;
use warnings;

use Genome;

# Notes from Chris Miller:
# bsub -oo err.log -q long
# -M 16000000 -R 'select[type==LINUX64 && mem>16000] span[hosts=1] rusage[mem=16000]' -n 8
# -J chimera -oo outputdir/chimera.err
# "python /gsc/bin/chimerascan_run.py -v -p 8
#   /gscmnt/sata921/info/medseq/cmiller/annotations/chimeraScanIndex/
#   $fastq1 $fastq2 $outputdir"

class Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult {
    is => "Genome::Model::RnaSeq::DetectFusionsResult",
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_output_directory();
    $self->_prepare_staging_directory();

    my $cmd_path = $self->_path_for_version($self->version);
    unless ($cmd_path) {
        die $self->error_message("Failed to find a path for chimerascan for version " . $self->version . "!");
    }

    my $index_dir = $self->_resolve_index_dir;

    my ($fastq1, $fastq2) = $self->_get_fastq_files_for_model();
    my $params = $self->detector_params;
    my $output_directory = $self->output_dir;

    my $n_threads = $self->_available_cpu_count;

    my $cmd = "python $cmd_path/chimerascan_run.py -v -p $n_threads $params $index_dir $fastq1 $fastq2 $output_directory >$output_directory/chimera_result.out";

    local $ENV{PYTHONPATH} =  ($ENV{PYTHONPATH} ? $ENV{PYTHONPATH} . ":" : "")  . $self->_python_path_for_version($self->version);

    $self->_put_bowtie_version_in_path();

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fastq1, $fastq2, $index_dir, $output_directory],
    );

    $DB::single=1;

    $self->_promote_data();
    $self->_remove_staging_directory();
    $self->_reallocate_disk_allocation();

    return $self;
}

sub _staging_disk_usage {
    #return enough to cover our temp dir (which will be under staging)
    return 60 * 1024 * 1024;
}

sub _path_for_version {
    my ($class,$version) = @_;
    die("You requested an unavailable version of Chimerascan. Requested: $version") unless $version eq '0.4.3';
    return $ENV{GENOME_SW} . "/chimerascan/chimerascan-$version/chimerascan";
}

sub _python_path_for_version {
    my ($class,$version) = @_;
    die("You requested an unavailable version of Chimerascan. Requested: $version") unless $version eq '0.4.3';
    return $ENV{GENOME_SW} . "/chimerascan/chimerascan-$version/build/lib.linux-x86_64-2.6";
}

sub _resolve_index_dir {
    my $self = shift;

    my $index = Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult::Index->get_or_create(
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        version => $self->version,
        bowtie_version => $self->alignment_result->bowtie_version,
        reference_build => $self->alignment_result->reference_build,
        annotation_build => $self->annotation_build,
    );

    unless($index){
        die("Unable to get a chimerascan index result");
    }

    return $index->output_dir;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/chimerascan/' . $self->id;
}

1;
