package Genome::Model::Tools::Soap::DaccDownload;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require Genome::Model::Tools::Dacc::Download;

class Genome::Model::Tools::Soap::DaccDownload {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
        version => {
            is => 'Text',
            doc => 'Version of Soap DeNovo to use',
            valid_values => [qw/ dacc /],
        },
        import_location => {
            is => 'Text',
            doc => 'The location of the assembly to import.',
        },
        output_dir_and_file_prefix => {
            is => 'Text',
            doc => 'Path and common prefix name for output files',
        },
    ],
};

sub help_brief {
    return 'Download Dacc soap denovo assembly';
}

sub help_detail {
    return <<HELP;
    Import a soap denovo assembly from the DACC.
HELP
}

sub execute {
    my $self = shift;

    $self->debug_message('Import SOAP assembly from the DACC');

    my $output_dir_and_file_prefix = $self->output_dir_and_file_prefix;
    $self->debug_message('Output directory and file prefix: '.$output_dir_and_file_prefix);
    my ($file_prefix, $output_dir) = File::Basename::fileparse($output_dir_and_file_prefix);
    if ( not -d $output_dir ) {
        $self->error_message("Invalid output directory: $output_dir does not exist");
        return;
    }
    my $edit_dir = $output_dir.'/edit_dir';
    mkdir $edit_dir if not -d $edit_dir;
    $self->debug_message('Output directory: '.$output_dir);
    $self->debug_message('Edit directory: '.$edit_dir);
    $self->debug_message('File prefix: '.$file_prefix);

    my @center_names = (qw/ Baylor LANL /);
    my ($center_name) = grep { $file_prefix =~ /$_/ } @center_names;
    if ( not $center_name ) {
        $self->error_message("Cannot determine center name from file prefix: $file_prefix");
        return;
    }
    $self->debug_message('Center name: '.$center_name);

    my $dacc_directory = $self->import_location;
    $self->debug_message('DACC directory: '.$dacc_directory);
    my $dacc_downloader = Genome::Model::Tools::Dacc::Download->create(
        dacc_directory => $dacc_directory,
        destination => $edit_dir,
    );
    if ( not $dacc_downloader ) {
        $self->error_message('Cannot creat DACC downloader.');
        return;
    }
    $dacc_downloader->dump_status_messages(1);

    my %available_files_and_sizes = $dacc_downloader->available_files_and_sizes;
    if ( not %available_files_and_sizes ) {
        $self->error_message('No files found in DACC directory: '.$dacc_directory);
        return;
    }

    $self->debug_message('Determining files to download');
    my $exts = {
        scafSeq => {
            required => 0,
            destination_dir => $output_dir,
        },
        agp=> {
            required => 1,
            destination_dir => $edit_dir,
        },
        'contigs.fa' => {
            required => 1,
            destination_dir => $edit_dir,
        },
        'scaffolds.fa' => {
            required => 1,
            destination_dir => $edit_dir,
        },
    };
    my @exts = keys %$exts;
    for my $ext ( @exts ) {
        my @available_files = grep { m/\.$ext/ } keys %available_files_and_sizes;
        if ( not @available_files ) {
            if ( $exts->{$ext}->{required} ) {
                $self->error_message("File ext ($ext) is required for download, but was not found in the DACC directory");
                return;
            }
            else {
                next;
            }
        }
        if ( @available_files == 1 ) {
            $exts->{$ext}->{file} = $available_files[0];
            next;
        }
        my ($available_pga_file) = grep { m/PGA/ } @available_files;
        if ( not $available_pga_file ) {
            $self->error_message("Found multiple files in dacc directory ($dacc_directory) for file ext ($ext), but one of them does not have PGA in the name. Files: @available_files");
            return;
        }
        $exts->{$ext}->{file} = $available_pga_file;
    }
    my @files_to_download = map { $exts->{$_}->{file} } grep { $exts->{$_}->{file} } @exts;
    $self->debug_message("Files to download: @files_to_download");

    $self->debug_message("Executing downloader");
    $dacc_downloader->files(\@files_to_download);
    if ( not $dacc_downloader->execute ) {
        $self->error_message('DACC downloader failed to execute');
        return;
    }
    $self->debug_message("Executing downloader...OK");

    $self->debug_message("Rename files");
    for my $ext ( @exts ) {
        next if not defined $exts->{$ext}->{file}; # OK
        my $file_name = $exts->{$ext}->{file};
        my $destination_dir = $exts->{$ext}->{destination_dir};
        if ( $file_name !~ /PGA/ and $destination_dir eq $edit_dir ) {
            next;
        }
        my $from = $edit_dir.'/'.$file_name;
        my $size = -s $from;
        my $to_base_name = $file_name;
        $to_base_name =~ s/PGA/$center_name/;
        my $to = $destination_dir.'/'.$to_base_name;
        $self->debug_message("Move $ext file");
        $self->debug_message("From $from to $to");
        $self->debug_message("Size: $size");
        my $move = File::Copy::move($from, $to);
        if ( not $move ) {
            $self->error_message('Move failed: '.$!);
            return;
        }
        my $new_size = -s $to;
        if ( not defined $new_size or $new_size != $size ) {
            $self->error_message("Move succeeded, but file ($to) now has different size: $size <=> ".(defined $new_size ? $new_size : 'undef'));
            return;
        }
        $self->debug_message("Move $ext file...OK");
    }
    $self->debug_message("Rename files...OK");

    $self->debug_message('Import SOAP assembly...OK');

    return 1;
}

1;

