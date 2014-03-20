package Genome::Model::GenotypeMicroarray::Build::CreateFilteredGenotypeTsvFile;

use strict;
use warnings;

use Genome;

use File::Basename;

class Genome::Model::GenotypeMicroarray::Build::CreateFilteredGenotypeTsvFile {
    is => 'Command::V2',
    has => {
        build => {
            is => 'Genome::Model::Build::GenotypeMicroarray',
            doc => 'Build to operate on.',
        },
    },
};

sub execute {
    my $self = shift; 
    $self->debug_message('Create filtered genotype TSV file...');
    $self->debug_message('Headers: not printed');
    $self->debug_message('Fields: chromosome, position, alleles');

    my $build = $self->build;
    my $genotype_file = $build->genotype_file_path;
    $self->debug_message('Filtered genotype file: '.$genotype_file);
    my $extract = Genome::Model::GenotypeMicroarray::Command::ExtractToCsv->create(
        build => $build,
        output => $genotype_file,
        fields => [qw/ chromosome position alleles /],
        headers => 0,
        filters => $build->model->genotype_filters,
    );
    if ( not $extract ) {
        $self->error_message('Failed to create command to create genotype file!');
        return;
    }

    $extract->dump_status_messages(1);

    if ( not $extract->execute ) {
        $self->error_message('Failed to execute command to create genotype file!');
        return;
    }

    if ( not -s $genotype_file ) {
        $self->error_message('Executed command to create genotype file, but file is empty! '.$genotype_file);
        return;
    }

    # Nutter made this file name, so we will link to it
    $self->debug_message('Link genotype file to gold2geno file...');
    $self->debug_message('Genotype file: '.$genotype_file);
    my $gold2geno_file = $build->gold2geno_file_path;
    $self->debug_message('Gold2geno file: '.$gold2geno_file);

    # Make a relative symlink if they are in the same directory. I think this
    # will always be the case but since gold2geno_file_path is not locally
    # defined I will check. Relative is better in case build's allocation is
    # moved or archived -> unarchived.
    if (dirname($genotype_file) eq dirname($gold2geno_file)) {
        Genome::Sys->create_symlink(basename($genotype_file), $gold2geno_file);
    } else {
        Genome::Sys->create_symlink($genotype_file, $gold2geno_file);
    }

    if ( not -l $gold2geno_file  or not -s $gold2geno_file ) {
        $self->error_message('Failed to link genotype file to gold2geno file!');
        return;
    }
    $self->debug_message('Link genotype file to gold2geno file...OK');

    $self->debug_message('Create filtered genotype TSV file...OK');
    return 1;
}

1;

