package Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard;

use strict;
use warnings;

use Genome;
use POSIX;

class Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard {
    is => 'Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler',
    has_optional => [
        max_jvm_heap_size => {
            is => 'Integer',
            doc => 'Size (in GB) of the JVM heap for Picard',
            is_constant => 1,
        },
    ],
};

sub max_gb { 32 };

sub default_max_jvm_heap_size {
    my $max_gb = max_gb();
    my $max_kb = 1_048_576 * $max_gb;
    my $default_max_jvm_heap_size = $max_gb;
    my $mem_limit_kb = Genome::Sys->mem_limit_kb;
    if ($mem_limit_kb) {
        my $safe_mem_limit_kb = int(0.8 * $mem_limit_kb);
        if ($max_kb > $safe_mem_limit_kb) {
            my $safe_mem_limit_gb = floor($safe_mem_limit_kb / 1_048_576);
            if ($safe_mem_limit_gb == 0) {
                die "Does not work on systems with less than 1GB of memory.\n";
            }
            $default_max_jvm_heap_size = $safe_mem_limit_gb;
        }
    }
    return $default_max_jvm_heap_size;
}

sub execute {
    my $self = shift;

    $self->max_jvm_heap_size($self->default_max_jvm_heap_size) unless $self->max_jvm_heap_size;

    my %mark_duplicates_params = (
        file_to_mark => $self->input_bam,
        marked_file => $self->output_path,
        metrics_file => $self->metrics_file,
        remove_duplicates => 0,
        tmp_dir => $self->scratch_directory,
        log_file => $self->log_file, 
        dedup_version => $self->version,
        dedup_params  => $self->parameters,
        max_jvm_heap_size => $self->max_jvm_heap_size,
    );

    if($self->include_comment) {
        $mark_duplicates_params{include_comment} = $self->include_comment;
    }

    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    my $mark_dup_cmd = Genome::Model::Tools::Sam::MarkDuplicates->create(%mark_duplicates_params);
    my $mark_dup_rv  = $mark_dup_cmd->execute;

    if ($mark_dup_rv != 1)  {
        $self->error_message("Error Marking Duplicates!");
        $self->error_message("Return value: ".$mark_dup_rv);
        $self->error_message("Check parameters and permissions in the RUN command above.");
        return;
    }

    return 1;
}

1;
