package Genome::InstrumentData::Command::Dacc::Status;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use File::Basename;
use File::Path;
use File::Copy;

class Genome::InstrumentData::Command::Dacc::Status {
    is  => 'Genome::InstrumentData::Command::Dacc',
};

sub help_brief {
    return 'Give the import status of a sample from the DACC';
}

sub help_detail {
    return help_brief();
}

sub execute {
    # Always return true from execute b/c this is just checking status of the Sra sample id
    my $self = shift;

    $self->debug_message('Status for '.$self->__display_name__);

    my $sample = $self->_get_sample;
    if ( not $sample ) { # no error in get
        # Not started
        $self->debug_message('Status: No Sample, needs download');
        return 1;
    }

    my @instrument_data = $self->_get_instrument_data;
    if ( not @instrument_data ) {
        $self->debug_message('Status: No instrument data, needs download');
        return 1;
    }

    my $have_archve_path = 0;
    for my $instrument_data ( @instrument_data ) {
        my $archive_path = eval{ $instrument_data->archive_path; };
        $have_archve_path++ if -e $archive_path;
    }

    if ( @instrument_data == $have_archve_path ) {
        $self->_update_library;
        $self->debug_message('Status: Done');
        return 1;
    }

    if ( not -d $self->_dl_directory ) {
        # Started, but failed b4 dl
        $self->debug_message('Status: No download directory, needs download');
        return 1;
    }

    my @existing_data_files = $self->existing_data_files;
    if ( not @existing_data_files ) {
        $self->debug_message('Status: No data files, needs download');
        return 1;
    }

    my $md5_ok = $self->_validate_md5;
    if ( not $md5_ok ) {
        $self->debug_message('Status: Failed to validate md5, needs download');
        return 1;
    }

    my $update_library = $self->_update_library;
    if ( not $update_library ) {
        $self->debug_message('Status: Failed to update library, xmls are corrupt');
        return 1;
    }

    # Everything seems ok...
    $self->debug_message('Status: Needs import');

    return 1
}

1;

