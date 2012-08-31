package Genome::InstrumentData::Command::Dacc::Import::Fastq;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::InstrumentData::Command::Dacc::Import::Fastq {
    is  => 'Genome::InstrumentData::Command::Dacc::Import',
    has => [
        format => {
            is => 'Text',
            is_constant => 1,
            value => 'fastq',
        },
        _paired_read_count => { is => 'Integer', is_optional => 1},
        _singleton_instrument_data => { is_optional => 1},
        _singleton_read_count => { is => 'Integer', is_optional => 1},
    ],
};

sub help_brief {
    return 'Import the downloaded fastqs from the DACC';
}

sub help_detail {
    return help_brief();
}

sub _downloaded_data_files {
    my $self = shift;

    my $dl_directory = $self->_dl_directory;
    my @downloaded_data_files = $self->existing_data_files;
    if ( not @downloaded_data_files ) {
        $self->error_message('No existing data files in download directory: '.$dl_directory);
        return;
    }

    my %files;
    my ($dl_fwd_fastq) = grep { m#\.1\.fastq# } @downloaded_data_files;
    if ( not -e $dl_fwd_fastq ) {
        $self->error_message('Forward fastq (*.1.fastq) not found in download directory: '.$dl_directory);
        return;
    }
    $files{fwd} = { file => $dl_fwd_fastq, new_name => 's_1_1_sequence.txt' }; 

    my ($dl_rev_fastq) = grep { m#\.2\.fastq# } @downloaded_data_files;
    if ( not -e $dl_rev_fastq ) {
        $self->error_message('Reverse fastq (*.2.fastq) not found in download directory: '.$dl_directory);
        return;
    }
    $files{rev} = { file => $dl_rev_fastq, new_name => 's_1_2_sequence.txt' }; 

    my ($dl_singleton_fastq) = grep { m#\.singleton\.fastq# } @downloaded_data_files;
    if ( not -e $dl_singleton_fastq ) {
        $self->error_message('Singleton fastq (*.singleton.fastq) not found in download directory: '.$dl_directory);
        return;
    }
    $files{singleton} = { file => $dl_singleton_fastq, new_name => 's_2_sequence.txt' }; 

    return \%files;
}

sub _execute {
    my $self = shift;

    my $unzip = $self->_unzip_fastqs;
    return if not $unzip;

    my $read_cnts = $self->_validate_fastq_read_counts;
    return if not $read_cnts;

    my $singleton = $self->_get_or_create_singleton_instrument_data;
    return if not $singleton;

    my $import = $self->_archive_and_update;
    return if not $import;

    return 1;
}

sub _unzip_fastqs {
    my $self = shift;

    my $dl_directory = $self->_dl_directory;
    my @zipped_fastqs = glob($dl_directory.'/*.fastq.bz2');

    return 1 if not @zipped_fastqs; # ok

    $self->status_message('Unzip fastqs...');

    for my $zipped_fastq ( @zipped_fastqs ) {
        my $cmd = "bunzip2 -f $zipped_fastq";
        $self->status_message($cmd);
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message("Cannot unzip fastq: $zipped_fastq");
            return;
        }
    }
    $self->status_message('Unzip fastqs...OK');

    return 1;
}

sub _validate_fastq_read_counts {
    my $self = shift;

    $self->status_message('Validate fastq read counts...');

    my $files = $self->_downloaded_data_files;
    my $fwd_read_count = $self->_read_count_for_fastq( $files->{fwd}->{file} );
    return if not defined $fwd_read_count;
    my $rev_read_count = $self->_read_count_for_fastq( $files->{rev}->{file} );
    return if not defined $rev_read_count;
    if ( $fwd_read_count != $rev_read_count ) {
        $self->error_message("Read counts for foward/reverse fastqs does not match: $fwd_read_count <=> $rev_read_count");
        return;
    }
    $self->_paired_read_count( $fwd_read_count + $rev_read_count );
    $self->status_message('Paired read count: '.$self->_paired_read_count);

    my $singleton_read_count = $self->_read_count_for_fastq( $files->{singleton}->{file} );
    return if not defined $singleton_read_count;
    $self->_singleton_read_count($singleton_read_count);
    $self->status_message('Singleton read count: '.$self->_singleton_read_count);

    $self->status_message('Validate fastq read counts...OK');

    return 1;
}

sub _read_count_for_fastq {
    my ($self, $fastq) = @_;

    my $line_count = `wc -l < $fastq`;
    if ( $? or not $line_count ) {
        $self->error_message("Line count on fastq ($fastq) failed : $?");
        return;
    }

    chomp $line_count;
    if ( ($line_count % 4) != 0 ) {
        $self->error_message("Line count ($line_count) on fastq ($fastq) not divisble by 4.");
        return;
    }

    return $line_count / 4;
}

sub _get_or_create_singleton_instrument_data {
    my $self = shift;

    my @instrument_data = $self->_get_instrument_data;
    if ( @instrument_data == 2 ) {
        $self->status_message('Got singelton instrument data: '.$instrument_data[1]->id);
        $self->_singleton_instrument_data($instrument_data[1]);
        return 1;
    }
    elsif ( @instrument_data > 2 ) {
        $self->error_message('Somehow there are more than two instrument data to import fastqs. Please fix.');
        return;
    }

    $self->status_message('Create singelton instrument data');

    my $singleton_fastq = $self->_singleton_fastq;
    my $size = -s $singleton_fastq;
    my $kilobytes_requested = int($size / 950); # 5% xtra space
    my $singleton_instrument_data = $self->_create_instrument_data(
        kilobytes_requested => $kilobytes_requested,
    );
    if ( not $singleton_instrument_data ) {
        $self->error_message('Failed to create singelton instrument data. See above errors.');
        return;
    }
    $self->_singleton_instrument_data($singleton_instrument_data);

    $self->status_message('Create singelton instrument data: '.$self->_singleton_instrument_data->id);

    return $singleton_instrument_data;
}

sub _archive_and_update {
    my $self = shift;

    $self->status_message('Archive and update...');

    my $paired_instrument_data = $self->_main_instrument_data;
    my $singleton_instrument_data = $self->_singleton_instrument_data;
    my $dl_directory = $self->_dl_directory;
    my $files = $self->_downloaded_data_files;

    #<PAIRED>#
    my $fwd_new_name = $files->{fwd}->{new_name};
    my $fwd_new_fastq = $dl_directory.'/'.$fwd_new_name;
    my $move = $self->_move_file($files->{fwd}->{file}, $fwd_new_fastq);
    return if not $move;
    my $rev_new_name = $files->{rev}->{new_name};
    my $rev_new_fastq = $dl_directory.'/'.$rev_new_name;
    $move = $self->_move_file($files->{rev}->{file}, $rev_new_fastq);
    return if not $move;
    my $archive = $self->_create_archive($paired_instrument_data, $fwd_new_name, $rev_new_name);
    return if not $archive;

    #<SINGLETON>#
    my $singleton_new_name = $files->{singleton}->{new_name};
    my $singleton_new_fastq = $dl_directory.'/'.$singleton_new_name;
    $move = $self->_move_file($files->{singleton}->{file}, $singleton_new_fastq);
    return if not $move;
    $archive = $self->_create_archive($singleton_instrument_data, $singleton_new_name);
    return if not $archive;

    $self->status_message('Update instrument data...');
    $paired_instrument_data->description($self->sra_sample_id.' Pairs that have been human screened, Q2 trimmed and deduped from the DACC');
    $paired_instrument_data->original_data_path(join(',', map { $files->{$_}->{file} } (qw/ fwd rev /)));
    $paired_instrument_data->read_count( $self->_paired_read_count );
    $paired_instrument_data->is_paired_end(1);

    $singleton_instrument_data->description($self->sra_sample_id.' Singletons that have been human screened, Q2 trimmed and deduped from the DACC');
    $singleton_instrument_data->original_data_path( $files->{singleton}->{file} );
    $singleton_instrument_data->read_count( $self->_singleton_read_count );
    $self->status_message('Update instrument data...OK');

    $self->status_message('Archive and update...OK');
    
    return 1;
}

sub _create_archive {
    my ($self, $instrument_data, @files) = @_;

    $self->status_message('Create archive for '.$instrument_data->id);

    Carp::confess('Cannot create archive b/c no instrument data was given') if not $instrument_data;
    Carp::confess('Cannot create archive b/c no files were given') if not @files;
    my $dl_directory = $self->_dl_directory;
    my @non_existing_files = grep { !-e $dl_directory.'/'.$_ } @files;
    Carp::confess('Cannot create archive b/c no these files do not exist: '.join(' ', @non_existing_files)) if @non_existing_files;

    my $temp_tar_file = $dl_directory.'/temp.tgz';
    unlink $temp_tar_file if -e $temp_tar_file;
    $self->status_message("Tar-ing files to $temp_tar_file");
    my $tar_cmd = "tar cvzf $temp_tar_file -C $dl_directory ".join(' ', @files);
    $self->status_message($tar_cmd);
    my $rv = eval { Genome::Sys->shellcmd(cmd => $tar_cmd); };
    if ( not $rv ) {
        $self->error_message("Tar command failed: $tar_cmd");
        return;
    }
    $self->status_message("Tar-ing OK");

    my $archive_path = $instrument_data->archive_path;
    unlink $archive_path if -e $archive_path;
    my $move = $self->_move_file($temp_tar_file, $archive_path);
    return if not $move;

    $self->status_message('Remove original files...');
    for my $file ( @files ) {
        unlink $file;
    }
    $self->status_message('Remove original files...OK');

    $self->status_message('Archive fastqs...OK');

    return 1;
}

1;

