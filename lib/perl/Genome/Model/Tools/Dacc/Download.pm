package Genome::Model::Tools::Dacc::Download;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;

class Genome::Model::Tools::Dacc::Download { 
    is => 'Genome::Model::Tools::Dacc',
    has => [
        destination => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'Directory to put downloaded files.',
        },
        files => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 3,
            doc => 'File names to download, space separated.'
        },
        launch_to_lsf => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Download must be run through LSF. This will schedule an LSF job with the correct rusage.',
        },
    ],
};

sub help_brief {
    return 'Download files from the DACC';
}

sub help_detail {
    return <<HELP;
    Download files from the DACC site. Give the DACC directory, files and the destination. 
    
    Because of the internet download resource requirement, this command must run through the LSF queue. It cannot be run locally. Specify --launch-to-lsf.
HELP
}

sub rusage {
    return (qw/ internet_download_mbps=100 /);
}

sub execute {
    my $self = shift;

    $self->status_message('Download');

    my $dacc_remote_directory = $self->dacc_remote_directory;
    $self->status_message("Dacc remote directory: $dacc_remote_directory");

    my @files = $self->files;
    my $files_string = join(' ', @files);
    $self->status_message("Files: $files_string");

    #This commented out because validation required sshing into the server
    #to determine that file exists but we can no longer ssh into ther server 5/13/11

    #my $files_exist = $self->validate_files_exist_in_dacc_directory(@files);
    #if ( not $files_exist ) {
    #    $self->error_message('Some/All file(s) for download do not exists on the dacc directory. See above errors.');
    #    return;
    #}

    my $destination = $self->destination;
    if ( not -d $destination ) {
        $self->error_message("Destination ($destination) does not exist or is not a direcory.");
        return;
    }
    $self->status_message("Destination: $destination");

    if ( $self->launch_to_lsf ) {
        return $self->_launch_to_lsf($destination, @files);
    }

    my $in_lsf = $self->validate_running_in_lsf_and_on_a_blade;
    return if not $in_lsf;

    my $remote_files_string = join(' ', map { $dacc_remote_directory.$_ } @files);
    my $cmd = $self->base_command.' '.$remote_files_string.' '.$destination;
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message("Aspera command failed: $cmd");
        return;
    }

    #This commented out because validation required sshing into the server
    #to determine that file exists but we can no longer ssh into ther server 5/13/11
    
    #my $dl_ok = $self->validate_files_were_downloaded;
    #return if not $dl_ok;

    $self->status_message('Download...OK');

    return 1;
}

sub validate_files_exist_in_dacc_directory {
    my $self = shift;

    $self->status_message('Validate files exist on the DACC');

    my %available_files_and_sizes = $self->available_files_and_sizes;
    if ( not %available_files_and_sizes ) {
        $self->error_message('No files availbale in the dacc directory: '. $self->dacc_directory);
        return;
    }

    my @files = $self->files;
    my $error;
    for my $file ( @files ) {
        if ( not exists $available_files_and_sizes{$file} ) {
            $error = 1;
            $self->error_message("File ($file) does not exist in the dacc directory: ". $self->dacc_directory);
        }
    }

    return if $error;

    $self->status_message('Validate files exist on the DACC...OK');

    return 1;
}

sub validate_files_were_downloaded {
    my $self = shift;

    $self->status_message('Validate files were downloaded');

    my %available_files_and_sizes = $self->available_files_and_sizes;
    if ( not %available_files_and_sizes ) {
        $self->error_message('No files availbale in the dacc directory: '. $self->dacc_directory);
        return;
    }

    my $destination = $self->destination;
    my @files = $self->files;
    my $error;
    for my $file_name ( @files ) {
        my $file = $destination.'/'.$file_name;
        my $size = -s $file;
        $self->status_message('File: '.$file);
        $self->status_message('Size: '.$size);
        $self->status_message('Size on DACC: '.$available_files_and_sizes{$file_name});
        if ( not $size ) {
            $error = 1;
            $self->error_message("Attempted to download file ($file_name) from the DACC, but it does not exist");
        }
        if ( not $size or $size != $available_files_and_sizes{$file_name} ) {
            $error = 1;
            $self->error_message("Attempted to download file ($file) from the DACC, but size is different: $size <=> $available_files_and_sizes{$file_name}");
        }
    }

    return if $error;

    $self->status_message('Validate files were downloaded...OK');

    return 1;
}

1;

=pod

=head1 Disclaimer

Copyright (C) 2005 - 2010 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

