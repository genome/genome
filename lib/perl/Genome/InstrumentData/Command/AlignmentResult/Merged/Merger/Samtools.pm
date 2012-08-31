package Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Samtools;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Samtools {
    is => 'Genome::InstrumentData::Command::AlignmentResult::Merged::Merger',
};

sub execute {
    my $self = shift;

    my $merge_cmd = Genome::Model::Tools::Sam::Merge->create(
        files_to_merge => [$self->input_bams],
        merged_file => $self->output_path,
        is_sorted => 1,
        #software => $merge_software,
        merger_name => 'samtools',
        merger_version => $self->version,
        merger_params => $self->parameters,
    );

    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->status_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    my $merge_rv = $merge_cmd->execute();

    if ($merge_rv != 1)  {
        $self->error_message("Error merging: ".join("\n", @{ $self->input_bams }));
        $self->error_message("Output target: " . $self->output_path);
        $self->error_message("Using software: samtools");
        $self->error_message("Version: ". $self->version);
        $self->error_message("You may want to check permissions on the files you are trying to merge.");
        return;
    }

    return 1;
}

1;
