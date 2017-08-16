package Genome::Model::SingleSampleGenotype::Command::CreateVcfTar;

use strict;
use warnings;

use Genome;

use File::Spec;

class Genome::Model::SingleSampleGenotype::Command::CreateVcfTar {
    is        => 'Command::V2',
    doc       => 'Create tar(s) of all vcfs for the build(s). One tar will be created for each build.',
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

        my $manifest_path = $self->build_manifest_file($build);

        my @tar_cmd = ( "tar", "-rvf", $tar_fullpath );

        for my $r ( sort { &result_sorter($a,$b) } $build->haplotype_caller_result ) {
            my $filename = $r->_vcf_filename;
            my $dir      = $r->output_dir;
            Genome::Sys->shellcmd( cmd => [ @tar_cmd, "-C", $dir, $filename ] );
        }

        my (undef, $manifest_dir, $manifest_file) = File::Spec->splitpath( $manifest_path );
        Genome::Sys->shellcmd( cmd => [ @tar_cmd, "-C", $manifest_dir, $manifest_file ] );

        Genome::Sys->shellcmd( cmd => [ "gzip", $tar_fullpath ] );
    }
    return 1;
}

sub build_manifest_file {
    my $self = shift;
    my $build = shift;

    my ($fh, $manifest_file) = Genome::Sys->create_temp_file('manifest.md5');

    for my $r ( sort { &result_sorter($a,$b) } $build->haplotype_caller_result ) {
        my $vcf = $r->vcf_file;
        my $filename = $r->_vcf_filename;

        my $md5 = Genome::Sys->md5sum($vcf);
        $fh->say(join("\t", $filename, $md5));
    }

    $fh->close;

    return $manifest_file;
}

sub result_sorter {
    return (($_[0]->intervals)[0] cmp ($_[1]->intervals)[0]);
}

sub sub_command_category {
    return 'analyst tools';
}

1;
