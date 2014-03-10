package Genome::Model::Tools::Mummer::Nucmer;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Mummer::Nucmer {
    is => 'Genome::Model::Tools::Mummer',
    has => [
        prefix => {
            is => 'Text',
            doc => 'Path to and prefix of output file',
        },
        reference => {
            is => 'Text',
            doc => 'Reference sequence',
        },
        query => {
            is => 'Text',
            doc => 'Query sequence',
        },
    ],
};

sub help_brief {
    'Tool to run nucmer aligner';
}

sub help_detail {
    return <<EOS
gmt mummer nucmer --prefix OUT --reference MRSA.MLST.fasta --query contigs.bases
EOS
}

sub execute {
    my $self = shift;

    if ( not -s $self->reference ) {
        $self->error_message('Failed to find reference file or file is zero size: '.$self->reference);
        return;
    }

    if ( not -s $self->query ) {
        $self->error_message('Failed to find query file or file is zero size: '.$self->query);
        return;
    }

    my $output_dir = File::Basename::dirname($self->prefix);
    if ( $output_dir and not -d $output_dir ) {
        $self->error_message("Failed to find output dir: $output_dir");
        return;
    }

    my $version_nucmer = $self->path_for_version.'/nucmer';
    if ( not -s $version_nucmer ) {
        $self->error_message('Failed to find nucmer in version, '.$self->path_for_version);
        return;
    }

    my $cmd = $version_nucmer.' --prefix '.$self->prefix.' '.$self->reference.' '.$self->query;

    $self->debug_message("Running nucmer with command: $cmd");

    my %cmd_params = (
        cmd          => $cmd,
        #input_files  => [ $self->reference, $self->query ],
        output_files => [ $self->prefix.'.delta' ],
    );
    my $rv = eval{ Genome::Sys->shellcmd( %cmd_params ); };
    if ( not $rv ) {
        $self->error_message("Failed to run numcer with command: $cmd");
        return;
    }

    return 1;
}

1;
