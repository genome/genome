package Genome::Model::Tools::RegulomeDb::ModifyRoisBasedOnScore;

use strict;
use warnings;
use Genome;
use WWW::Mechanize;

class Genome::Model::Tools::RegulomeDb::ModifyRoisBasedOnScore {
    is => 'Genome::Model::Tools::RegulomeDb',
    has => [
        roi_list => {
            is => 'String',
            doc => 'Path to roi list in bed format (0-based)',
        },
        scored_regions => {
            is => 'String',
            doc => 'Bed file with scored regions (0-based)',
        },
        output_file => {
            is => 'String',
            doc => 'Path of output file',
        },
        valid_scores => {
            is => 'String',
            is_many => 1,
            doc => 'Filter ROIs, keeping only these scores',
            valid_values => [qw(1 2 3 4 5 6 1a 1b 1c 1d 1e 1f 2a 2b 2c 3a 3b)],
        },
    ],
};

sub execute {
    my $self = shift;
    #get only scored regions with good scores
    my $good_scores = Genome::Sys->create_temp_file_path;
    my $in = Genome::Sys->open_file_for_reading($self->scored_regions);
    my $out = Genome::Sys->open_file_for_writing($good_scores);
    while (my $line = <$in>) {
        chomp $line;
        if ($self->good_score($line)) {
            $out->print("$line\n");
        }
    }
    $out->close;

    #Get only regions of ROIs that overlap with good scored-regions
    Genome::Model::Tools::BedTools::Intersect->execute(
        input_file_a => $self->roi_list,
        input_file_b => $good_scores,
        intersection_type => 'overlaps',
        output_file => $self->output_file,
        input_file_a_format => 'bed',
    );

    return 1;
}

sub good_score {
    my $self = shift;
    my $line = shift;
    my @fields = split(/\t/, $line);
    return grep {$fields[3] =~ /$_/} $self->valid_scores;
}

1;

