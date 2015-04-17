package Genome::InstrumentData::Command::Import::WorkFlow::SourceFile;

use strict;
use warnings;

use Genome;

require File::Basename;
require List::MoreUtils;
require LWP::Simple;

class Genome::InstrumentData::Command::Import::WorkFlow::SourceFile { 
    is => 'UR::Object',
    has => {
        source_file => { is => 'Text', },
    },
    has_transient => {
        format => { is =>'Text', },
    },
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    die $self->error_message('No source file given!') if not $self->source_file;

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
    my $source_file = $self->source_file;
    $source_file =~ s/\Q.$_\E$// for Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs->types;
    my ($source_file_base_name, $path, $suffix) = File::Basename::fileparse(
        $source_file, keys %suffixes_to_original_format
    );
    if ( not $suffix or not $suffixes_to_original_format{$suffix} ) {
        die $self->error_message('Unrecognized source file format! '.$source_file_base_name);
    }

    my $format = $suffixes_to_original_format{$suffix};
    if ( $self->is_tar ) {
        die $self->error_message("Cannot process tar $format! %s", $source_file) if $format ne 'fastq';
        $format .= '_archive';
    }

    return $self->format($format);
}

sub retrieval_method {
    my $self = shift;

    if ( $self->source_file =~ m#^(http|https|file)\://# ) {
        return 'remote url';
    }
    else {
        return 'local disk';
    }
}

sub file_size {
    my $self = shift;

    my $source_file = $self->source_file;
    if ( $source_file !~ /^http/ ) {
        return -s $self->source_file;
    }

    my $agent = LWP::UserAgent->new;
    my $response = $agent->head($source_file);
    if ( not $response->is_success ) {
        $self->error_message($response->message) if $response->message;
        die $self->error_message('HEAD failed for remote file! %s', $source_file);
    }

    return $response->headers->content_length;
}

sub is_tar {
    my $self = shift;
    my $source_file = $self->source_file;
    return List::MoreUtils::any { $source_file =~ /\Q$_\E$/ } (qw/ .tgz .tar .tar.gz /);
}

sub kilobytes_required_for_processing {
    my $self = shift;

    my %formats_and_multipliers = (
        bam => 3,
        fastq => 2,
        fastq_archive => 2,
        sra => 4,
    );

    my $source_file = $self->source_file;
    my $size = $self->file_size($source_file);
    die $self->error_message('Source file does not have any size! '.$source_file) if not $size;

    my $kb_required = int($size / 1024) + 1; #convert to kb, +1 for ceiling
    $kb_required *= 3 if $source_file =~ /\.gz$/; # assume ~30% compression rate for gzipped fastq

    my $multiplier = $formats_and_multipliers{ $self->format };
    $kb_required *= $multiplier; # for conversion to bam

    return $kb_required;
}

1;

