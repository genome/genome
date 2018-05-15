package Genome::InstrumentData::Command::Import::Inputs::SourceFile;

use strict;
use warnings;

use Genome;

require File::Basename;
require List::MoreUtils;
require LWP::Simple;

class Genome::InstrumentData::Command::Import::Inputs::SourceFile { 
    is => 'UR::Object',
    has => {
        path => { is => 'Text', },
    },
    has_transient => {
        format => { is =>'Text', },
    },
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    die $self->error_message('No source file given!') if not $self->path;

    $self->_resolve_format;

    return $self;
}

sub _resolve_format {
    my $self = shift;

    my %suffixes_to_original_format = (
        fastq => 'fastq',
        fq => 'fastq',
        bam => 'bam',
        sra => 'sra',
    );
    my $source_file = $self->path;
    $source_file =~ s/\.gz$//;
    my ($source_file_base_name, $path, $suffix) = File::Basename::fileparse(
        $source_file, keys %suffixes_to_original_format
    );
    if ( not $suffix or not $suffixes_to_original_format{$suffix} ) {
        die $self->error_message('Unrecognized source file format! '.$source_file_base_name);
    }

    return $self->format( $suffixes_to_original_format{$suffix} );
}

sub retrieval_method {
    my $self = shift;

    if ( $self->path =~ m#^(http|https|file)\://# ) {
        return 'remote url';
    }
    else {
        return 'local disk';
    }
}

sub _is_file_local {
    my $self = shift;

    if ( $self->path =~ m#^(http|https|file)\://# ) {
        return;
    }

    return 1;
}

sub file_size {
    return $_[0]->_file_size($_[0]->path);
}

sub _file_size {
    my ($self, $file) = @_;

    if ( $self->_is_file_local($file) ) {
        return -s $file;
    }

    my $agent = LWP::UserAgent->new;
    my $response = $agent->head($file);
    if ( not $response->is_success ) {
        $self->error_message($response->message) if $response->message;
        die $self->error_message('HEAD failed for remote file! %s', $file);
    }

    return $response->headers->content_length;
}

sub is_tar {
    my $self = shift;
    my $source_file = $self->path;
    return List::MoreUtils::any { $source_file =~ /\Q$_\E$/ } (qw/ .tgz .tar .tar.gz /);
}

sub kilobytes_required_for_processing {
    my $self = shift;

    my %formats_and_multipliers = (
        bam => 3,
        fastq => 2,
        sra => 4,
    );

    my $source_file = $self->path;
    my $size = $self->file_size($source_file);
    die $self->error_message('Source file does not have any size! '.$source_file) if not $size;

    my $kb_required = int($size / 1024) + 1; #convert to kb, +1 for ceiling
    $kb_required *= 3 if $source_file =~ /\.gz$/; # assume ~30% compression rate for gzipped fastq

    my $multiplier = $formats_and_multipliers{ $self->format };
    $kb_required *= $multiplier; # for conversion to bam

    return $kb_required;
}

sub md5_path {
    return $_[0]->path.'.md5';
}

sub md5_path_size {
    return $_[0]->_file_size($_[0]->md5_path);
}

1;

