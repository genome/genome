package Genome::Model::Tools::RefCov::Topography;

use strict;
use warnings;

use Genome;
use GD::Graph::lines;

class Genome::Model::Tools::RefCov::Topography {
    is => ['Command'],
    has => [
        bam_file => {
            is => 'Text',
        },
        bed_file => {
            is => 'Text',
        },
        min_depth_filter => {
            is_optional => 1,
            default_value => 1,
        },
        wingspan => {
            is_optional => 1,
            default_value => 0,
        },
        graph => {
            is => 'Boolean',
            doc => 'Creates a png format line graph',
            default_value => 0,
        },
        output_dir => {
            is_optional => 1,
            default_value => '.',
        }
    ],
};

sub help_brief {
    "Generate a topography of coverage or a graph using GD.",
}

sub help_synopsis {
my $self = shift;
    return <<"EOS"
gmt ref-cov topography...
EOS

};

sub help_detail {
    '
Given a BAM file and a list of ROI in a BED file, this command will generate
either a text topography of coverage depth per position or a GD graph of
the topography of each ROI.
'
}

sub execute {
    my $self = shift;

    my $alignments  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->bam_file );
    unless ($alignments) {
        die('Failed to load alignment file '. $self->bam_file);
    }
    my $regions = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $self->bed_file,
        wingspan => $self->wingspan,
    );
    unless ($regions) {
        die('Failed to load BED file '. $self->bed_file );
    }
    my $bam = $alignments->bio_db_bam;
    my $index = $alignments->bio_db_index;
    my $output_dir = $self->output_dir;
    while (my $region = $regions->next_region) {
        my $tid = $alignments->tid_for_chr($region->{chrom});
        my $coverage = $index->coverage( $bam, $tid, ($region->{start} - 1), $region->{end});
        unless (scalar( @{ $coverage } ) == $region->{length}) {
            die('The length of region '. $region->{name} .' '. $region->{id}
                    .'('. $region->{length} .') does not match the coverage array length '. scalar( @{ $coverage }));
        }
        my $topology = Genome::Model::Tools::RefCov::Topology->create(
            coverage => $coverage,
            min_depth => $self->min_depth_filter,
        );
        if ($self->graph) {
            my $output_file = $output_dir .'/'. $region->{name} .'.png';
            my $graph = GD::Graph::lines->new();
            $graph->set(
                'x_label' => 'Base Pair Position',
                'x_label_skip' => 1000,
                'y_label' => 'Read Depth',
                'title' => $region->{name} .' Topology',
                marker_size => 1,
            );
            my $start = $region->{start};
            my $stop = $region->{end};
            my @positions = ($start .. $stop);
            my @data = (\@positions,$topology->topology);
            my $gd = $graph->plot(\@data);
            open(IMG, '>'. $output_file) or die $!;
            binmode IMG;
            print IMG $gd->png;
            close IMG;
        } else {
            $topology->print_topology;
        }
    }
    return 1;
}

1;
