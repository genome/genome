package Genome::Model::RnaSeq::DetectFusionsResult::DeFuseResult;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::DeFuseResult{
    is => "Genome::Model::RnaSeq::DetectFusionsResult",
    has_param => [
        detector_params => {
            doc => 'params to pass along to deFuse including (required) ensembl_release=num'
        }
    ]
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_output_directory();
    $self->_prepare_staging_directory();

    my $config_params = $self->detector_params;
    $config_params =~ s/ensemble_release=(\d+)\s*/ /i;
    my $ensembl_release = $1;

    die($self->error_message("You must specify an ensembl_release version in the detector params")) unless $ensembl_release;

    my $index_result = Genome::Model::RnaSeq::DetectFusionsResult::DeFuseResult::Index->get_or_create(
        ensembl_release => $ensembl_release,
        config_params => $config_params,
        version => $self->version
    );

    #regenerate the same config file as before, but use the software results staging dir
    $index_result->generate_config_file($self->temp_staging_directory);
    $self->run_defuse($index_result);

    $self->_promote_data();
    $self->_remove_staging_directory();
    $self->_reallocate_disk_allocation();

    return $self;
}

sub run_defuse {
    my $self = shift;
    my $index_result = shift;

    my $cmd = Genome::Model::RnaSeq::DetectFusionsResult->_path_for_command($self->version, "defuse.pl");
    $cmd .= " -c " . $self->temp_staging_directory . "/config.txt";
    $cmd .= " -d " . $index_result->output_dir;
    $cmd .= " -o " . $self->temp_staging_directory;
    $cmd .= " -p " . $self->_available_cpu_count;

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->temp_staging_directory . "/config.txt"],
        output_files => [map {$self->temp_staging_directory . "/$_" }
            qw{ results.txt results.filtered.txt results.classify.txt } ]
    );
}

sub _staging_disk_usage {
    return 60 * 1024 * 1024;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/de-fuse/' . $self->id;
}

1;
