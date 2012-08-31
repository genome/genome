package Genome::InstrumentData::Command::Dacc::Download;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require File::Path;

class Genome::InstrumentData::Command::Dacc::Download {
    is  => 'Genome::InstrumentData::Command::Dacc',
    has => [
        import_log_file => {
            is => 'Text',
            is_optional => 1,
            #shell_args_position => 3,
            doc => 'Log file to use when launching import job. If none given will use email.',
        },
    ],
};

#< Helps >#
sub help_brief {
    return 'Download and import from the DACC';
}

sub help_detail {
    return help_brief();
}
#<>#

#< Execute >#
sub execute {
    my $self = shift;

    $self->status_message('Download '.$self->sra_sample_id.' '.$self->format);

    if ( not $self->_is_host_a_blade ) {
        $self->error_message('To download from the DACC, this command must be run on a blade');
        return;
    }

    my $sample = $self->_get_sample;
    if ( not $sample ) {
        $sample = $self->_create_sample;
    }
    return if not $sample;

    my @instrument_data = $self->_get_instrument_data;
    if ( not @instrument_data ) {
        my $new_instrument_data = $self->_create_instrument_data(
            kilobytes_requested => $self->kb_to_request,
        );
        return if not $new_instrument_data;
        @instrument_data = $self->_get_instrument_data;
        return if not @instrument_data;
    }

    my @instrument_data_w_archives = grep { -e } map { eval{ $_->archive_path} } @instrument_data;
    if ( @instrument_data_w_archives ) {
        $self->error_message('Sra sample ('.$self->__display_name__.') has already been imported. These instrument data have a data file: '.join(' ', map { $_->id } @instrument_data_w_archives));
        return;
    }

    my $download = $self->_run_aspera;
    return if not $download;

    my $dl_dir_ok = $self->_dl_directory_exists;
    return if not $dl_dir_ok;

    my $md5_ok = $self->_validate_md5;
    return if not $md5_ok;

    $self->_launch_import; # no error check
    $self->_main_instrument_data->description('downloaded');

    $self->status_message('Download...OK');

    return 1;
}

sub _is_host_a_blade {
    my $self = shift;

    my $hostname = `hostname`;
    if ( not defined $hostname ) {
        $self->error_message('Cannot get hostname');
        return;
    }
    
    return $hostname =~ /blade/ ? 1 : 0;
}

sub _run_aspera {
    my $self = shift;

    $self->status_message('Run aspera...');

    my $dl_directory = $self->_dl_directory;
    if ( -d $dl_directory ) {
        # we tried dl'ing already - blow this dir away, ignore errors
        $self->status_message('Remove previous download directory: '.$dl_directory);
        File::Path::rmtree($dl_directory);
        Carp::confess('Tried to remove download directory, but could not') if -d $dl_directory;
        $self->status_message('Removing previous download directory...OK');
    }

    my $absolute_path = $self->_absolute_path;
    Carp::confess("Main absolute path does not exist: $absolute_path") if not -d $absolute_path;
    my $user = 'jmartin';
    my $sra_sample_id = $self->sra_sample_id;
    my $key_file = '/gsc/scripts/share/certs/dacc/dacc.ppk';
    if ( not -s $key_file ) {
        $self->error_message("Aspera key file ($key_file) does not exist.");
        return;
    }
    my $dacc_location = $self->dacc_location;

    # Exclude other formats 
    my $exclude = join(' ', map { '-E *'.$_.'*' } grep { $self->format ne $_ } $self->valid_formats);

    my $cmd = "ascp -QTd -l100M -i $key_file $exclude $user\@aspera.hmpdacc.org:$dacc_location/$sra_sample_id ".$absolute_path;
    $self->status_message($cmd);
    my $rv = eval { Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message("Aspera command failed: $cmd");
        return;
    }

    $self->status_message('Run aapera...OK');

    return 1;
}
#<>#

#< Launch Import >#
sub _launch_import {
    my $self = shift;

    my $sra_sample_id = $self->sra_sample_id;
    my $sub_command_format = join('-', split(' ', $self->format));
    my $logging;
    if ( defined $self->import_log_file ) {
        unlink $self->import_log_file if -e $self->import_log_file;
        $logging = '-oo '.$self->import_log_file;
    }
    else {
        $logging = '-u '. Genome::Config->user_email;
    }

    my $cmd = "bsub -q long $logging genome instrument-data dacc import $sub_command_format $sra_sample_id";
    if ( not $self->validate_md5 ) {
        $cmd .= ' --novalidate-md5';
    }
    $self->status_message('Launch import: '.$cmd);
    my $rv = eval { Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message('Failed to launch import: '.$@);
    }

    return 1;
}

1;

