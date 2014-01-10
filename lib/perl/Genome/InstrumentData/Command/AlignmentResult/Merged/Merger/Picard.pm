package Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard;

use strict;
use warnings;

use Genome;
use POSIX;

class Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard {
    is => 'Genome::InstrumentData::Command::AlignmentResult::Merged::Merger',
    has_optional => [
        max_jvm_heap_size => {
            is => 'Integer',
            doc => 'Size (in GB) of the JVM heap for Picard',
            is_constant => 1,
        },
    ],
};

sub max_gb { 12 };

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

    my $merge_cmd = Genome::Model::Tools::Sam::Merge->create(
        files_to_merge => [$self->input_bams],
        merged_file => $self->output_path,
        is_sorted => 1,
        bam_index => 0,
        merger_name => 'picard',
        merger_version => $self->version,
        merger_params  => $self->parameters,
        use_version => $self->samtools_version,
        max_jvm_heap_size => $self->max_jvm_heap_size,
        include_comment => $self->include_comment,
    );

    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    my $merge_rv = $merge_cmd->execute();

    if ( not $merge_rv )  {
        $self->error_message("Error merging: ".join("\n", $self->input_bams));
        $self->error_message("Output target: " . $self->output_path);
        $self->error_message("Using software: picard");
        $self->error_message("Version: ". $self->version);
        $self->error_message("You may want to check permissions on the files you are trying to merge.");
        return;
    }

    return 1;
}

1;
