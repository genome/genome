package Genome::Model::SingleSampleGenotype::Command::CreateVcfTar;

use strict;
use warnings;

use Genome;

class Genome::Model::SingleSampleGenotype::Command::CreateVcfTar {
    is        => 'Command::V2',
    doc       => 'Create a tar of all vcfs for the build',
    has_input => {
        builds => {
            is      => "Genome::Model::Build::SingleSampleGenotype",
            is_many => 1,
        },
        output_dir => {
            is  => 'DirectoryPath',
            doc => 'The directory for the created tars.'
        }
    },
};

sub execute {
    my $self = shift;

    for my $build ( $self->builds ) {
        my $output_dir   = $self->output_dir;
        my $build_id     = $build->build_id;
        my $tar_fullpath = File::Spec->join( $output_dir, "$build_id.tar" );

        if ( -e $tar_fullpath ) {
            $self->fatal_message("Tar file already exists, please remove: $tar_fullpath");
        }

        my @tar_cmd = ( "tar", "-rvf", $tar_fullpath );

        for my $r ( sort { ( $a->intervals )[0] cmp( $b->intervals )[0] } $build->haplotype_caller_result ) {
            my $filename = $r->_vcf_filename;
            my $dir      = $r->output_dir;
            Genome::Sys->shellcmd( cmd => [ @tar_cmd, "-C", $dir, $filename ] );
        }

        Genome::Sys->shellcmd( cmd => [ "gzip", $tar_fullpath ] );
    }
    return 1;
}

sub sub_command_category {
    return 'analyst tools';
}

1;
