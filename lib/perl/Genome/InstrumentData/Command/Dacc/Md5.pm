package Genome::InstrumentData::Command::Dacc::Md5;

use strict;
use warnings;

use Genome;

require Cwd;
use Data::Dumper 'Dumper';

class Genome::InstrumentData::Command::Dacc::Md5 {
    is => 'Command',
    has => [
        data_files => {
            is => 'Text',
            is_many => 1,
            doc => 'Data files to validate with MD5',
        },
        format => {
            is => 'Text',
            doc => 'Format of the SRA data to validate.',
        },
        confirmed_md5_file => {
            is => 'Text',
            doc => 'MD5 file for confirmation. If does not exist, MD5s will be generated and saved.',
        },
        # priv
        _sra_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'SRA sample id',
        },
        _directory => {
            is => 'Text',
            is_optional => 1,
            doc => 'Directory where the data and md5 files live.',
        },
        _data_file_names => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'Data file base names',
        },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    if ( not $self->format ) {
        $self->error_message('No format given to validate MD5.');
        return;
    }

    if ( not $self->data_files ) {
        $self->error_message('No data files given to validate MD5.');
        return;
    }

    my $confirmed_md5_file = $self->confirmed_md5_file;
    if ( not $confirmed_md5_file ) {
        $self->error_message('No confirmed Md5 file given');
        return;
    }

    my %directories;
    my @data_file_names;
    for my $data_file ( $self->data_files ) {
        my ($name, $dir) = File::Basename::fileparse($data_file);
        $dir =~ s/\/$//;
        $directories{$dir} = 1;
        push @data_file_names, $name;
    }

    my @directories = keys %directories;
    if ( @directories > 1 ) {
        $self->error_message('Data files are in multiple directories: '.Dumper(\@directories));
        $self->delete;
        return;
    }
    $self->_directory($directories[0]);
    $self->_data_file_names(\@data_file_names);
    $self->_sra_id( File::Basename::basename($self->_directory) );

    return $self;
}

sub execute {
    return $_[0]->_validate;
}

# Files
sub expected_md5_files { 
    my $self = shift;
    my $method = '_expected_md5_files_for_'.$self->format;
    return $self->$method;
}

sub _expected_md5_files_for_fastq {
    my $self = shift;
    return ( $self->_sra_id.'.md5' => [ map { s/\.bz2$//; $_ } $self->_data_file_names ] );
}

sub _expected_md5_files_for_bam {
    my $self = shift;
    my %expected_md5s;
    for my $data_file_name ( $self->_data_file_names ) {
        $expected_md5s{ $data_file_name.'.md5' } = $data_file_name;
    }
    return %expected_md5s;
}

sub _expected_md5_files_for_sff {
    my $self = shift;
    my %expected_md5s;
    for my $data_file_name ( $self->_data_file_names ) {
        my $md5_file_name = $data_file_name;
        $md5_file_name =~ s/sff$/md5/;
        $expected_md5s{$md5_file_name} = $data_file_name;
    }
    return %expected_md5s;
}

# Validate
sub _validate {
    my $self = shift;

    $self->status_message('Validate md5...');

    my %dacc_md5 = $self->_load_dacc_md5;
    return if not %dacc_md5;

    my %confirmed_md5 =$self->_load_confirmed_md5;
    return if not %confirmed_md5;

    for my $file ( keys %confirmed_md5 ) {
        if ( not defined $dacc_md5{$file} ) {
            $self->error_message("No MD5 from the DACC for file ($file)");
            return;
        }
        if ( $dacc_md5{$file} ne $confirmed_md5{$file} ){
            $self->error_message("MD5 for file ($file) does not match: $dacc_md5{$file} <=> $confirmed_md5{$file}");
            return;
        }
    }

    $self->status_message('Validate md5...OK');

    return 1;
}

# Load
sub _load_dacc_md5 {
    my $self = shift;

    $self->status_message('Load DACC md5...');

    my $cwd = Cwd::cwd();
    my $directory = $self->_directory;
    chdir $directory or die "Cannot change dir to $directory";

    my %md5_files = $self->expected_md5_files;
    Carp::confess('No expected md5 files.') if not %md5_files;

    my %files_and_md5;
    for my $md5_file ( keys %md5_files ) {
        # md5_file => file(s)
        if ( not -e $md5_file ) {
            $self->error_message("Expected MD5 file ($md5_file) does not exist.");
            return;
        }
        my %current_files_and_md5 = $self->_load_md5($md5_file);
        return if not %current_files_and_md5;

        if ( ref $md5_files{$md5_file} ) {
            # one md5 file has many md5s - use the file name in the md5 file to associate
            # file names should match here, and exist
            for my $file ( @{$md5_files{$md5_file}} ) {
                $files_and_md5{$file} = $current_files_and_md5{$file};
            }
        }
        else {
            # the md5 for the file is the value in the current hash, but the file name may not be the same
            my $file = $md5_files{$md5_file}; 
            my @md5s = values %current_files_and_md5;
            if ( @md5s > 1 ) { # should only be one
                $self->error_message("Expected one md5 in file ($md5_file) to match $file, but got: ".Dumper(\%current_files_and_md5));
                return;
            }
            $files_and_md5{$file} = $md5s[0];
        }
    }

    chdir $cwd;

    $self->status_message('Load DACC md5...OK');

    return %files_and_md5;
}

sub _load_confirmed_md5 {
    my $self = shift;

    $self->status_message('Load confirmed md5...');

    my $md5_file = $self->confirmed_md5_file;
    if ( not -e $md5_file ) {
        my $generate = $self->_generate;
        return if not $generate;
    }

    my %files_and_md5 = $self->_load_md5($md5_file);
    if ( not %files_and_md5 ) {
        $self->error_message("Cannot load confirmed md5 from file: $md5_file");
        return;
    }

    $self->status_message('Load confirmed md5...OK');

    return %files_and_md5;
}

sub _load_md5 {
    my ($self, $md5_file) = @_;

    my $md5_fh = eval{ Genome::Sys->open_file_for_reading($md5_file); };
    if ( not defined $md5_fh ) {
        $self->error_message("Cannot open md5 file ($md5_file): $@");
        return;
    }

    my $sra_id = $self->_sra_id;
    my %files_and_md5;
    while ( my $line = $md5_fh->getline ) {
        chomp $line;
        my ($md5, $file) = split(/\s+/, $line);
        $file =~ s#^$sra_id/##;
        $file =~ s/\.bz2$//;
        $files_and_md5{$file} = $md5;
    }

    return %files_and_md5;
}

# Generate
sub _generate {
    my $self = shift;

    $self->status_message('Generate md5...');

    my $md5_file = $self->confirmed_md5_file;
    unlink $md5_file if -e $md5_file;

    my $cwd = Cwd::cwd();
    my $directory = $self->_directory;
    chdir $directory;

    for my $file ( sort { $a cmp $b } $self->data_files ) {
        my $file_name = File::Basename::basename($file);
        my $cmd = "md5sum $file_name >> $md5_file";
        $self->status_message("MD5 command: $cmd");
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message("Failed to run md5sum on $file: $@");
            return;
        }
    }

    chdir $cwd;

    $self->status_message('Generate md5...OK');

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

