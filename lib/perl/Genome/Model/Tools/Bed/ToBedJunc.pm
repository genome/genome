package Genome::Model::Tools::Bed::ToBedJunc;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::Bed::ToBedJunc {
    is => ['Command::V2'],
    has_input => [
        bed12_file => {
            is => 'Text',
            doc => 'The input BED12 format file of junctions.',
        },
        bed_file => {
            is => 'Text',
            doc => 'The output junctions file in BED format.',
        }
    ],
};

sub help_brief {
    "Tool to convert BED12 files to BED junction files.",
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

sub execute {
    my $self = shift;

    my $bed12_reader = Genome::Utility::IO::BedReader->create(
        input => $self->bed12_file,
        headers => Genome::Utility::IO::BedReader->bed12_headers,
    );
    my $bed_writer = Genome::Utility::IO::BedWriter->create(
        output => $self->bed_file,
    );
    while (my $bed12_data = $bed12_reader->next) {
        unless (
            ($bed12_data->{'chrom'} =~ /\w+/) &&
                ($bed12_data->{'chromStart'} =~ /\d+/) &&
                    ($bed12_data->{'chromEnd'} =~ /\d+/)
                ) {
            next;
        }

        my $blocks = $bed12_data->{'blockCount'};
        if ($blocks == 1) { next; }

        my $start = $bed12_data->{'chromStart'};
        my @b_sizes = split(/,/, $bed12_data->{'blockSizes'});
        my @b_offsets = split(/,/, $bed12_data->{'blockStarts'});
        for (my $b = 1; $b < $blocks; $b++) {
            my $left = $start + $b_offsets[$b-1] + $b_sizes[$b-1];
            my $right = $start + $b_offsets[$b] + 1;
            my %bed_data = (
                'chr' => $bed12_data->{chrom},
                'start' =>  $left,
                'end' => $right,
                'name' => $bed12_data->{name},
                'score' => $bed12_data->{score},
                'strand' => $bed12_data->{'strand'},
            );
            $bed_writer->write_one(\%bed_data);
        }
    }
    return 1;
}


1;

