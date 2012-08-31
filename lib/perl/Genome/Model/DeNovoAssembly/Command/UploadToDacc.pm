package Genome::Model::DeNovoAssembly::Command::UploadToDacc; 

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require File::Temp;
require IO::File;

class Genome::Model::DeNovoAssembly::Command::UploadToDacc {
    is => 'Genome::Command::Base',
    has => [
        model => {
            is => 'Genome::Model',
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'Filter text to get de novo assembly model(s).',
        },
        _sra_sample_id => { is => 'Text', is_optional => 1, },
    ],
};

#< Help >#
sub help_brief {
    return 'Upload to the DACC';
}

sub help_detail {
    return help_brief();
}
#<>#

#< Execute >#
sub execute {
    my $self = shift;

    my $build = $self->_get_last_succeeeded_build_from_model;
    return if not defined $build;

    my $verify = $self->_verify_build($build);
    return if not $verify;

    my $upload = $self->_upload_files($build);
    return if not $upload;

    return 1;
}
#<>#

#< Files Names Exts >#
sub _files_and_exts {
    return (
        soap_scaffold_sequence_file => 'scafSeq',
        pga_agp_file => 'agp',
        pga_contigs_fasta_file => 'contigs.fa',
        pga_scaffolds_fasta_file => 'scaffolds.fa',
    );
}

sub _file_names_to_upload {
    my %files_and_exts = _files_and_exts();
    return keys %files_and_exts;
}

sub _get_last_succeeeded_build_from_model {
    my $self = shift;

    my $model = $self->model;
    if ( not $model->isa('Genome::Model::DeNovoAssembly') ) {
        $self->error_message('Model ('.$model->id.') is not a de novo assembly model. It is "'.$model->type_name.'". Cannot upload to the DACC.');
        return;
    }

    if ( $model->processing_profile->assembler_base_name ne 'soap' ) {
        $self->error_message('Model ('.$model->id.') was not assembled with soap. It was assembled with '.$model->processing_profile->assembler_base_name.'. Cannot upload to the DACC.');
        return;
    }

    my $build = $model->last_succeeded_build;
    if ( not defined $build ) {
        $self->error_message('Model ('.$model->id.') does not have a last succeeded build. Cannot upload to the DACC.');
        return;
    }

    return $build;
}

sub _verify_build {
    my ($self, $build) = @_;

    $self->status_message('Verify assembly length...');
    my $assembly_length = $build->assembly_length;
    if ( not defined $assembly_length ) {
        $self->error_message('Assembly length is not defiend for last succeeded build: '.$build->id.'. Cannot upload to the DACC.');
        return;
    }
    if ( $assembly_length == 0 ) {
        $self->error_message('Assembly length is 0 for last succeeded build: '.$build->id.'. Cannot upload to the DACC.');
        return;
    }
    $self->status_message('Verify assembly length OK: '.$assembly_length);

    $self->status_message('Verify files to upload...');
    for my $file_name ( _file_names_to_upload() ) {
        my $file = $build->$file_name;
        if ( not -s $file ) {
            $self->error_message("File $file for $file_name does not have any size.");
            return;
        }
    }
    $self->status_message('Verify files to upload...OK');

    $self->status_message('Get sra sample id...');
    my %sra_sample_ids = map { $_->sra_sample_id => 1 } grep { defined $_->sra_sample_id } $build->instrument_data;
    if ( not %sra_sample_ids ) {
        $self->error_message('No sra sample ids found in instrument data for build: '.$build->id);
        return;
    }
    my @sra_sample_ids = keys %sra_sample_ids;
    if ( @sra_sample_ids > 1 ) {
        $self->error_message('Multiple sra sample ids found in build instrument data: '.$build->id);
        return;
    }
    $self->_sra_sample_id($sra_sample_ids[0]);
    $self->status_message('Got sra sample id: '.$self->_sra_sample_id);

    return 1;
}

sub _upload_files {
    my ($self, $build) = @_;

    # sra id
    my $sra_sample_id = $self->_sra_sample_id;
    die "No sra sample id" if not defined $sra_sample_id;

    # md5 file
    $self->status_message('Open md5 file...');
    my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
    my $md5_file = $tmpdir.'/'.$sra_sample_id.'.md5';
    my $md5_fh = IO::File->new($md5_file, 'w');
    if ( not defined $md5_fh ) {
        $self->error_message('Cannot open md5 file: '.$md5_file);
        return;
    }
    $self->status_message('Open md5 file...OK');

    # go through files, get md5, push to upload
    $self->status_message('Upload files...');
    my %files_and_exts = _files_and_exts();
    my @files_to_upload = ( $md5_file );
    for my $file_name ( sort keys %files_and_exts ) {
        my $file = $build->$file_name;
        $self->status_message("Checking $file_name...");
        if ( not -s $file ) {
            $self->error_message("File $file for $file_name does not have any size.");
            return;
        }
        push @files_to_upload, $file;
        $self->status_message("Checking $file_name...OK: $file");

        $self->status_message('MD5 for '.$file);
        my $md5 = $self->_md5_for_file($file);
        return if not defined $md5;
        my $base_name = File::Basename::basename($file);
        my $md5_string = $md5."\t".$sra_sample_id.'_WUGC/'.$base_name;
        $self->status_message($md5_string);
        $md5_fh->print("$md5_string\n");
        $self->status_message('MD5 OK: '.$md5);
    }
    $md5_fh->flush;
    my $aspera = $self->_run_aspera(@files_to_upload);
    return if not $aspera;
    $self->status_message('Upload files...OK');

    return 1;
}

sub _run_aspera {
    my ($self, @files) = @_;

    my $sra_sample_id = $self->_sra_sample_id;
    die "No sra sample id" if not defined $sra_sample_id;
    my $dest_dir = '/WholeMetagenomic/03-Assembly/PGA/'.$sra_sample_id.'_WUGC/';

    $self->status_message('Aspera key file...');
    my $key_file = '/gscuser/ebelter/DACC/keys/dacc.ppk';
    if ( not -s $key_file ) {
        $self->error_message("Asepeara key file ($key_file) deos not exist.");
        return;
    }
    $self->status_message('Aspera key file...OK');

    my $files = join(' ', @files);
    my $cmd = "ascp -QTd -l100M -i $key_file $files jmartin\@aspera.hmpdacc.org:$dest_dir";
    $self->status_message("Aspera: $cmd");
    #my $rv = 1;
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message("Aspera command failed: $cmd");
        return;
    }
    $self->status_message("Aspera OK");

    return 1;
}

sub _md5_for_file {
    my ($self, $file) = @_;

    if ( not -e $file ) {
        $self->error_message("Cannot get md5 for non existing file: $file");
        return;
    }

    my ($md5) = split(/\s+/, `md5sum $file`);

    if ( not defined $md5 ) {
        $self->error_message("No md5 returned for file: $file");
        return;
    }

    return $md5;
}

1;

#$HeadURL$
#$Id$
