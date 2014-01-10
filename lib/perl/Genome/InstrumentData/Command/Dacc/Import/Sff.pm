package Genome::InstrumentData::Command::Dacc::Import::Sff;

use strict;
use warnings;

use Genome;

require Carp;
use Data::Dumper 'Dumper';
require File::Copy;

class Genome::InstrumentData::Command::Dacc::Import::Sff {
    is  => 'Genome::InstrumentData::Command::Dacc::Import',
    has => [
        format => {
            is => 'Text',
            is_constant => 1,
            value => 'sff',
        },
    ],
};

sub help_brief {
    return 'Import the downloaded SFFs from the DACC';
}

sub help_detail {
    return help_brief();
}

sub _execute {
    my $self = shift;

    my @sffs = $self->existing_data_files;
    $self->debug_message('SFF count: '.scalar(@sffs));
    my @instrument_data;
    # Get the inst data that need the archive path, as the command may have died
    for my $instrument_data ( $self->_get_instrument_data ) {
        my $archive_path = eval{ $instrument_data->archive_path; };
        next if $archive_path and -e $archive_path;
        push @instrument_data, $instrument_data
    }

    if ( @instrument_data > @sffs ) {
        $self->error_message('Somehow there are more instrument data w/o archive paths ('.scalar(@instrument_data).') than SFFs ('.scalar(@sffs).') to move. Please fix.');
        return;
    }

    for ( my $i = 0; $i <= $#sffs; $i++ ) {
        my $instrument_data = $instrument_data[$i];
        my $sff = $sffs[$i];
        my $size = -s $sff;
        my $kilobytes_requested = int($size / 950); # 5% xtra space
        # Make sure we got an inst data for this sff
        if ( not $instrument_data ) { 
            $instrument_data = $self->_create_instrument_data(kilobytes_requested => $kilobytes_requested);
            if ( not $instrument_data ) {
                $self->error_message('Cannot create instrument data.');
                return;
            }
        }
        # Make sure the inst data has an allocation
        my $allocation = $instrument_data->allocations;
        if ( not $allocation ) { # this could happen if a command gets killed
            $allocation = $self->_create_instrument_data_allocation(
                instrument_data => $instrument_data,
                kilobytes_requested => $kilobytes_requested
            );
            return if not $allocation;
        }

        # move sff to archive path
        my $archive_path = eval{ $instrument_data->archive_path; };
        if ( not $archive_path ) {
            $self->error_message('No archive path for instrument data: '.$instrument_data->id);
            return;
        }
        $self->debug_message("Move $sff to $archive_path");
        my $move_ok = File::Copy::move($sff, $archive_path);
        if ( not $move_ok ) {
            $self->error_message("Failed to move SFF $sff to $archive_path: $!");
            return;
        }
        if ( not -e $archive_path ) {
            $self->error_message('Move succeeded, but archive path does not exist.');
            return;
        }
        if ( $size != -s $archive_path ) {
            $self->error_message("Moved SFF $sff to $archive_path but now file size is different.");
            return;
        }

        # update inst data, commit
        $instrument_data->original_data_path($sff);
        my $sff_base_name = File::Basename::basename($sff);
        $instrument_data->description($self->sra_sample_id." SFF $sff_base_name from the DACC");
        if ( not UR::Context->commit ) {
            $self->error_message('Cannot commit after moving SFF and uipdating instrument data.');
            return;
        }
    }

    return 1;
}

1;

