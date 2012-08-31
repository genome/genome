package Genome::Model::Tools::Dacc::Upload;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::Dacc::Upload { 
    is => 'Genome::Model::Tools::Dacc',
    has => [
        files => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 2,
            doc => '',
        },
        launch_to_lsf => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Download must be run through LSF. This will schedule an LSF job with the correct rusage.',
        },
    ],
};

sub help_brief {
    return 'Upload files to the DACC';
}

sub help_detail {
    return <<HELP;
    Upload files to the DACC site. Give the DACC directory and files. Once uploaded, it will be verified that they exist in the DACC directory, and the size is the same. 
    
    Because of the internet upload and aspera upload resource requirement, this command must run through the LSF queue. It cannot be run locally. Specify --launch-to-lsf.
HELP
}

sub rusage {
    return (qw/ internet_upload_mbps=100 aspera_upload_mbps=100 /);
}

sub execute {
    my $self = shift;

    $self->status_message('Upload');

    my $dacc_remote_directory = $self->dacc_remote_directory;
    $self->status_message("Dacc remote directory: $dacc_remote_directory");

    my @files = $self->files;
    for my $file ( @files ) {
        if ( not -e $file ) {
            $self->error_message("File to upload ($file) does not exist");
            return;
        }
    }
    my $files_string = join(' ', @files);
    $self->status_message("Files: $files_string");

    if ( $self->launch_to_lsf ) {
        return $self->_launch_to_lsf(@files);
    }

    my $in_lsf = $self->validate_running_in_lsf_and_on_a_blade;
    return if not $in_lsf;

    my $cmd = $self->base_command.' -d '.$files_string.' '.$dacc_remote_directory;
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message("Aspera command failed: $cmd");
        return;
    }
    $self->status_message('Aspera command OK');

    #This commented out because validation required sshing into the server
    #to determine that file exists but we can no longer ssh into ther server 5/13/11
    #my $upload_ok = $self->validate_files_were_uploaded;
    #return if not $upload_ok;

    $self->status_message('Upload...OK');

    return 1;
}

sub validate_files_were_uploaded {
    my $self = shift;

    $self->status_message('Validate files were uploaded');

    my $dacc_directory = $self->dacc_directory;
    my @files = $self->files;
    my %available_files_and_sizes = $self->available_files_and_sizes;

    if ( not %available_files_and_sizes ) {
        $self->error_message("No files found on dacc directory($dacc_directory). Odd, since the aspera command successfully completed.");
        return
    }

    my $error;
    for my $file ( @files ) {
        $self->status_message('File: '.$file);
        my $size = -s $file;
        $self->status_message('Size: '.$size);
        my $file_name = File::Basename::basename($file);
        Carp::confess("Cannot get base name for file: $file") if not $file_name;
        $self->status_message('File name: '.$file_name);
        if ( not exists $available_files_and_sizes{$file_name} ) {
            $error = 1;
            $self->error_message("Attempted to upload file ($file), but it is not in the DACC directory: $dacc_directory");
            next;
        }
        $self->status_message('Size on DACC: '.$available_files_and_sizes{$file_name});
        if ( $size != $available_files_and_sizes{$file_name} ) {
            $error = 1;
            $self->error_message("Attempted to upload file ($file), but file size in DACC directory ($dacc_directory) does not match: $size <=> $available_files_and_sizes{$file_name}");
        }
    }

    return if $error;

    $self->status_message('Validate files were uploaded...OK');

    return 1;
}

1;

=pod

=head1 Disclaimer

Copyright (C) 2005 - 2009 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

