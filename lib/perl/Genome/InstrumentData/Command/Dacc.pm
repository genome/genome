package Genome::InstrumentData::Command::Dacc;

use strict;
use warnings;

use Genome;

require File::Copy;
use Data::Dumper 'Dumper';
require IPC::Run;

class Genome::InstrumentData::Command::Dacc {
    is  => 'Command',
    is_abstract => 1,
    has => [
        sra_sample_id => {
            is => 'Text',
            is_input => 1,
            shell_args_position => 1,
            doc => 'SRA id to download and import from the DACC.',
        },
        format => {
            is => 'Text',
            is_input => 1,
            shell_args_position => 2,
            valid_values => [ valid_formats() ],
            doc => 'Format of the SRA id to download.',
        },
        validate_md5 => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Validate MD5 for data files.',
        },
        # sample
        _sample => { is_optional => 1, },
        _library => { is_optional => 1, },
        # inst data
        _main_instrument_data => { is_optional => 1, },
        _allocation => { via => '_main_instrument_data', to => 'disk_allocations' },
        _absolute_path => { via => '_allocation', to => 'absolute_path' },
        # dl directory
        _dl_directory => {
            calculate_from => [qw/ _absolute_path sra_sample_id /], 
            calculate => q| return $_absolute_path.'/'.$sra_sample_id |,
        },
        _dl_directory_exists => {
            calculate_from => [qw/ _dl_directory /], 
            calculate => q| return -d $_dl_directory ? $_dl_directory : undef; |,
        },
    ],
};

#< COMMAND >#
sub help_brief {
    return 'Donwload and import from the DACC';
}

sub help_detail {
    return help_brief();
}

sub __display_name__ {
    return $_[0]->sra_sample_id.' '.$_[0]->format;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    if ( $self->sra_sample_id !~ /^SRS/ ) {
        $self->error_message('Invalid sra id: '. $self->sra_sample_id);
        return;
    }

    return $self;
}
#<>#

#< FORMATS >#
sub formats_and_info {
    return (
        fastq => {
            sequencing_platform => 'solexa',
            import_format => 'sanger fastq',
            dacc_location => '/WholeMetagenomic/02-ScreenedReads/ProcessedForAssembly',
            kb_to_request => 100_000_000, # 100 GiB
        },
        sff => {
            sequencing_platform => '454',
            dacc_location => '/WholeMetagenomic/02-ScreenedReads/ProcessedForAssembly',
            kb_to_request => 15_000_000, # 15 GiB 
        },
        bam => {
            sequencing_platform => 'solexa',
            dacc_location => '/WholeMetagenomic/05-Analysis/ReadMappingToReference',
            kb_to_request => 20_000_000, # 20 GiB
        }
    );
}

sub valid_formats {
    my %formats = formats_and_info();
    return sort { $a cmp $b} keys %formats;
}

sub import_format {
    my $self = shift;
    my %formats = formats_and_info();
    return $formats{ $self->format }->{import_format} || $self->format;
}

sub sequencing_platform {
    my $self = shift;
    my %formats = formats_and_info();
    return $formats{ $self->format }->{sequencing_platform};
}

sub dacc_location {
    my $self = shift;
    my %formats = formats_and_info();
    return $formats{ $self->format }->{dacc_location};
}

sub kb_to_request {
    my $self = shift;
    my %formats = formats_and_info();
    return $formats{ $self->format }->{kb_to_request};
}

sub existing_data_files {
    my $self = shift;

    my $dl_directory = $self->_dl_directory_exists;
    return if not $dl_directory;

    my $format = $self->format;
    return sort { $a cmp $b } grep { $_ !~ /md5|xml/i } glob($dl_directory.'/*'.$self->format.'*');
}
#<>#

#< SAMPLE >#
sub _get_sample {
    my $self = shift;

    $self->debug_message('Get sample...');

    my $sample = Genome::Sample->get(name => $self->sra_sample_id);
    return if not defined $sample;
    $self->_sample($sample);
    $self->debug_message('Sample: '.join(' ',  map { $sample->$_ } (qw/ id name /)));

    my $library = $self->_get_or_create_library;
    return if not $library;
    return if not $self->_library;

    return $self->_sample;
}

sub _create_sample {
    my $self = shift;

    my $sra_sample_id = $self->sra_sample_id;
    if ( $sra_sample_id !~ /^SRS/ ) {
        $self->error_message("Invalid sra sample id: $sra_sample_id");
        return;
    }

    $self->debug_message("Create sample for $sra_sample_id");

    my $sample = Genome::Sample->create(
        name => $sra_sample_id,
        extraction_label => $self->sra_sample_id,
        cell_type => 'unknown',
    );

    if ( not defined $sample ) {
        $self->error_message("Cannot create sample for $sra_sample_id");
        return;
    }
    if ( not UR::Context->commit ) {
        $self->error_message('Cannot commit sample');
        return;
    }
    $self->_sample( $sample );
    $self->debug_message('Sample: '.join(' ',  map { $sample->$_ } (qw/ id name /)));

    my $library = $self->_get_or_create_library;
    return if not $library;

    $self->debug_message('Create sample...OK');

    return $self->_sample;
}

sub _get_or_create_library {
    my $self = shift;

    my $sample = $self->_sample;
    my $library_name = $sample->name.'-extlibs';
    my $library = Genome::Library->get(name => $library_name);
    if ( not $library ) {
        $library = Genome::Library->create(
            name => $library_name,
            sample_id => $sample->id,
        );
        if ( not $library ) {
            $self->error_message('Cannot create library: '.$library_name);
            return;
        }
        if ( not UR::Context->commit ) {
            $self->error_message('Cannot commit library: '.$library_name);
            return;
        }
    }
    $self->_library($library);

    $self->debug_message('Library: '.join(' ',  map { $self->_library->$_ } (qw/ id name /)));

    return $self->_library;
}
#<>#

#< INST DATA >#
sub _get_instrument_data {
    my $self = shift;

    if ( not $self->_library ) {
        Carp::confess('No library set! Cannot get instrument data w/o it!');
    }

    my @instrument_data = Genome::InstrumentData::Imported->get(
        library_id => $self->_library->id,
        import_source_name => 'DACC',
        import_format => $self->import_format,
    );
    return if not @instrument_data;
    $self->_main_instrument_data($instrument_data[0]);

    if ( not $self->_allocation ) { # main allocation for downloading
        my $allocation = $self->_create_instrument_data_allocation(
            instrument_data => $instrument_data[0],
            kilobytes_requested => $self->kb_to_request,
        );
        return if not $allocation;
    }

    $self->debug_message('Instrument data: '.join(' ', map { $_->id } @instrument_data));
    $self->debug_message('Absolute path: '.$self->_absolute_path );

    return @instrument_data;
}

sub _create_instrument_data {
    my ($self, %params) = @_;

    Carp::confess('No kb to request when creating instrument data.') if not $params{kilobytes_requested};

    $self->debug_message('Create instrument data...');

    my @instrument_data = $self->_get_instrument_data;
    my $instrument_data = Genome::InstrumentData::Imported->create(
        library_id => $self->_library->id,
        sra_sample_id => $self->sra_sample_id, 
        sequencing_platform => $self->sequencing_platform,
        import_format => $self->import_format,
        import_source_name => 'DACC',
        original_data_path => 0, 
        description => 'new',
        subset_name => scalar(@instrument_data) + 1,
    );
    if ( not $instrument_data ) {
        $self->error_message('Cannot create instrument data for sra sample id: '.$self->sra_sample_id);
        return;
    }
    if ( not UR::Context->commit ) {
        $self->error_message('Cannot commit instrument data.');
        return;
    }

    $params{instrument_data} = $instrument_data;
    my $allocation = $self->_create_instrument_data_allocation(%params);
    return if not $allocation;

    $self->debug_message('Create instrument data: '.$instrument_data->id);

    return $instrument_data;
}

sub _create_instrument_data_allocation {
    my ($self, %params) = @_;

    $self->debug_message('Create instrument data allocation...');

    my $instrument_data = $params{instrument_data};
    Carp::confess('No instrument data given to create allocation') if not $instrument_data;

    my $kilobytes_requested = $params{kilobytes_requested};
    Carp::confess('No kilobytes requested given to create allocation') if not $kilobytes_requested;

    my $allocation = $instrument_data->allocations;
    if ( defined $allocation ) {
        $self->debug_message('Allocation already exists for instrument data: '.$instrument_data->id);
        return;
    }

    $allocation = Genome::Disk::Allocation->allocate(
        owner_id => $instrument_data->id,
        owner_class_name => $instrument_data->class,
        disk_group_name => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
        allocation_path => 'instrument_data/imported/'.$instrument_data->id,
        kilobytes_requested => $kilobytes_requested,
    );

    if ( not $allocation ) {
        $self->error_message('Could not create disk allocation for instrument data: '.$instrument_data->id);
        return;
    }

    if ( not $allocation ) {
        $self->error_message('No allocation for instrument data: '.$instrument_data->id);
    }

    $self->debug_message('Create instrument data allocation...OK');

    return $allocation;
}
#<>#

#< MD5 >#
sub _validate_md5 {
    my $self = shift;

    $self->debug_message('Validate MD5...');

    if ( not $self->validate_md5 ) {
        $self->debug_message('Skip validate md5...');
        return 1;
    }

    my @data_files = $self->existing_data_files;
    if ( not @data_files ) {
        $self->error_message('No data files found for format: '.$self->format);
        return;
    }

    my $md5 = Genome::InstrumentData::Command::Dacc::Md5->create(
        data_files => \@data_files,
        format => $self->format,
        confirmed_md5_file => $self->_dl_directory.'/confirmed.md5',
    );
    if ( not $md5 ) {
        $self->error_message('Cannot create MD5 validate object.');
        return;
    }
    $md5->dump_status_messages(1);
    if ( not $md5->execute ) {
        $self->error_message('Failed to validate MD5');
        return;
    }

    return 1;
}
#<>#

#< Update Library >#
sub _update_library {
    my $self = shift;

    $self->debug_message('Update library...');

    my $dl_directory = $self->_dl_directory;
    if ( not -d $dl_directory ) {
        $self->error_message("Download directory ($dl_directory) does not exist.");
        return;
    }

    my @xml_files = glob($dl_directory.'/*.xml');
    if ( not @xml_files ) { # ok
        $self->debug_message('Attempt to update library, but no XMLs in download directory. This is OK.');
        return 1;
    }

    my $update_library = Genome::Sample::Command::Import::Dacc->create(
        sra_sample_id => $self->sra_sample_id,
        xml_files => \@xml_files,
    );
    if ( not $update_library ) {
        $self->error_message('Cannot create update library object.');
        return;
    }
    $update_library->dump_status_messages(1);
    if ( not $update_library->execute ) {
        $self->error_message('Failed to update library');
        return;
    }

    return 1;
}
#<>#

#< FILES ON THE DACC SITE >#
sub data_files_on_the_dacc {
    my $self = shift;

    my $oringinal_cert = $self->original_certificate;

    my ($in, $out);
    my $harness = IPC::Run::harness(['ssh', '-i', '/home/archive/dacc.sshkey', 'jmartin@aspera.hmpdacc.org',], \$in, \$out);
    $harness->pump until $out;

    $out = '';
    $in = "ls\n";
    $harness->pump until $out;
    print $out;

    return 1;
}
#<>#

#< MOVE FILE >#
sub _move_file {
    my ($self, $file, $new_file) = @_;

    Carp::confess("Cannot move file b/c none given!") if not $file;
    Carp::confess("Cannot move $file b/c it does not exist!") if not -e $file;
    Carp::confess("Cannot move $file to $new_file b/c no new file was given!") if not $new_file;

    $self->debug_message("Move $file to $new_file");

    my $size = -s $file;
    $self->debug_message("Size: $size");

    my $move_ok = File::Copy::move($file, $new_file);
    if ( not $move_ok ) {
        $self->error_message("Failed to move $file to $new_file: $!");
        return;
    }

    if ( not -e $new_file ) {
        $self->error_message('Move succeeded, but archive path does not exist.');
        return;
    }

    if ( $size != -s $new_file ) {
        $self->error_message("Moved $file to $new_file but now file size is different.");
        return;
    }

    $self->debug_message("Move...OK");

    return 1;
}
#<>#

1;

