package Genome::Model::GenotypeMicroarray::Build::CreateCopyNumberTsvFile;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Build::CreateCopyNumberTsvFile {
    is => 'Command::V2',
    has_input_output => {
        build => { is => 'Genome::Model::Build::GenotypeMicroarray', },
    },
};

sub execute {
    my $self = shift;
    $self->debug_message('Create copy number TSV file...');
    $self->debug_message('Headers: not printed');
    $self->debug_message('Fields: chromosome, position, log r ratio');

    my $build = $self->build;
    my $copy_number_file = $build->copy_number_file_path;
    $self->debug_message('Copy number file: '.$copy_number_file);
    my $extract = Genome::Model::GenotypeMicroarray::Command::ExtractToCsv->create(
        build => $build,
        output => $copy_number_file,
        fields => [qw/ chromosome position log_r_ratio /],
        headers => 0,
        filters => $build->model->genotype_filters,
    );
    if ( not $extract ) {
        $self->error_message('Failed to create command to create copy number file!');
        return;
    }

    $extract->dump_status_messages(1);

    if ( not $extract->execute ) {
        $self->error_message('Failed to execute command to create copy number file!');
        return;
    }
    if ( not -s $copy_number_file ) {
        $self->error_message('Executed command to create copy number file, but file is empty! '.$copy_number_file);
        return;
    }

    $self->debug_message('Create copy number file...OK');
    return 1;
}

1;

