package Genome::Model::Tools::Dacc::TarAndUpload;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require Genome::Model::Tools::Dacc::Upload;

class Genome::Model::Tools::Dacc::TarAndUpload { 
    is => 'Genome::Model::Tools::Dacc',
    has => [
        tar_file => {
            is => 'Bolean',
            shell_args_position => 3,
            doc => 'Tar file.',
        },
        files => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 4,
            doc => 'Files to tar and zip.',
        },
    ],
};

sub help_brief {
    return 'Tar files and upload to the DACC';
}

sub help_detail {
    return <<HELP;
    Tar files and then upload to the DACC site. Give the DACC directory and files. Once the files are tarred, the upload command will be executed. The original files and tar file will not be deleted.
HELP
}

sub execute {
    my $self = shift;

    my $tar_file = $self->tar_file;
    $self->status_message("Tar file: $tar_file");
    if ( $tar_file !~ /\.tgz$/ and $tar_file !~ /\.tar\.gz/  ) {
        $self->error_message('Tar file must have .tar.gz or .tgz extension');
        return;
    }

    my $file_string = join(' ', $self->files);
    $self->status_message("Files: $file_string");
    for my $file ( $self->files ) {
        if ( not -e $file ) {
            $self->error_message("File to upload ($file) does not exist");
            return;
        }
    }

    my $cmd = "tar cvzf $tar_file $file_string";
    my $rv = eval { Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message("Tar command failed: $cmd");
        return;
    }
    if ( not -e $tar_file ) {
        $self->error_message('Tar command succeeded, but no tar file was created');
        return;
    }
    $self->status_message("Tar-ing...OK");

    my $upload = Genome::Model::Tools::Dacc::Upload->create(
        dacc_directory => $self->dacc_directory,
        files => [ $tar_file ],
        launch_to_lsf => 1,
    );
    if ( not $upload ) {
        $self->error_message("Failed to launch upload, but tar file exists: $tar_file");
        return;
    }
    $upload->dump_status_messages(1);
    if ( not $upload->execute ) {
        $self->error_message("Failed to launch upload, but tar file exists: $tar_file");
        return;
    }

    return 1;
}

1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2009 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

