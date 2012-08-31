package Genome::Model::Tools::Bed::ToJunc;

use warnings;
use strict;

use Genome;

my @JUNC_HEADERS = qw/
                       chr:start-end
                       read_count
                     /;

class Genome::Model::Tools::Bed::ToJunc {
    is => ['Command::V2'],
    has_input => [
        bed_file => {
            is => 'Text',
            doc => 'The BED12 format file to convert to a junctions file.',
        },
        bed_format => {
            is => 'Text',
            is_optional => 1,
            doc => 'The BED format of the input file.',
            valid_values => ['bed12','bed'],
            default_value => 'bed',
        },
        junc_file => {
            is => 'Text',
            doc => 'The output junctions format file.',
        }
    ],
};

sub help_brief {
    "Tool to convert BED files to junction files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed to-junc...
EOS
}

sub help_detail {
    return <<EOS
    This script will transform tophat's junction.bed file to actual junction coordinates
EOS
}

sub default_junc_headers {
    return \@JUNC_HEADERS;
}

sub execute {
    my $self = shift;
    if ($self->bed_format eq 'bed12') {
        $self->convert_bed12;
    } elsif ($self->bed_format eq 'bed') {
        $self->convert_bed;
    }
    return 1;
}

sub convert_bed12 {
    my $self = shift;
    
    my $bed_reader = Genome::Utility::IO::BedReader->create(
        input => $self->bed_file,
        headers => Genome::Utility::IO::BedReader->bed12_headers,
    );
    unless ($bed_reader) {
        $self->error_message('Failed to load BED12 reader for file: '. $self->bed_file);
        die($self->error_message);
    }
    my $junc_writer = $self->_load_junc_writer;
    while (my $bed_data = $bed_reader->next) {
        unless (
            ($bed_data->{'chrom'} =~ /\w+/) &&
                ($bed_data->{'chromStart'} =~ /\d+/) &&
                    ($bed_data->{'chromEnd'} =~ /\d+/)
                ) {
            next;
        }

        my $blocks = $bed_data->{'blockCount'};
        if ($blocks == 1) { next; }

        my $start = $bed_data->{'chromStart'};
        my @b_sizes = split(/,/, $bed_data->{'blockSizes'});
        my @b_offsets = split(/,/, $bed_data->{'blockStarts'});
        for (my $b = 1; $b < $blocks; $b++) {
            my $left = $start + $b_offsets[$b-1] + $b_sizes[$b-1];
            my $right = $start + $b_offsets[$b] + 1;
            my $strand = $bed_data->{'strand'};
            my %junc_data = (
                'chr:start-end' => $bed_data->{chrom} .':'. $left .'-'. $right .'('. $strand .')',
                'read_count' => $bed_data->{score},
            );
            $junc_writer->write_one(\%junc_data);
        }
    }
}

sub convert_bed {
    my $self = shift;

    my $bed_reader = Genome::Utility::IO::BedReader->create(
        input => $self->bed_file,
    );
    unless ($bed_reader) {
        $self->error_message('Failed to load BED reader for file: '. $self->bed_file);
        die($self->error_message);
    }
    my $junc_writer = $self->_load_junc_writer;
    while (my $bed_data = $bed_reader->next) {
        my %junc_data = (
            'chr:start-end' => $bed_data->{chr} .':'. $bed_data->{start} .'-'. $bed_data->{end} .'('. $bed_data->{strand} .')',
            'read_count' => $bed_data->{score},
        );
        $junc_writer->write_one(\%junc_data);
    }
}

sub _load_junc_writer {
    my $self = shift;
    
    my $junc_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->junc_file,
        separator => "\t",
        headers => $self->default_junc_headers,
    );
    unless ($junc_writer) {
        $self->error_message('Failed to load junction writer for file: '. $self->junc_file);
        die($self->error_message);
    }
    return $junc_writer;
}
