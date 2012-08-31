package Genome::Model::Tools::Kmer::PlotDuplicationRate;

use strict;
use warnings;

use Genome;
use GD::Graph::lines;

class Genome::Model::Tools::Kmer::PlotDuplicationRate {
    is => ['Command'],
    has => [
        occratio_file => {
            doc => 'The occratio output file to plot',
        },
        plot_file => {
            doc => 'The path to write output png plot',
        },
    ],
};

sub execute {
    my $self = shift;

    my $occratio_fh = Genome::Sys->open_file_for_reading($self->occratio_file);
    my @x;
    my @y;
    while (my $line = $occratio_fh->getline){
        chomp($line);
        if ($line =~ /^#/) { next; }
        my ($mer_size,$count,$ratio) = split(/\s+/,$line);
        push @x, $mer_size;
        push @y, $ratio;
    }
    $occratio_fh->close;
    my $data = [\@x,\@y];
    my $graph = new GD::Graph::lines(800,600);
    $graph->set(
                x_label         => 'K-mer Size',
                y_label         => 'Duplication Rate',
                title           => 'K-mer Based Duplication Rate',
                y_max_value     => 1,
                x_label_position => .5,
                x_label_skip    => 5,
            )
        or warn $graph->error;
    my $gd = $graph->plot($data) or die $graph->error;
    open(IMG, '>'. $self->plot_file) or die $!;
    binmode IMG;
    print IMG $gd->png;
    close(IMG);
    return 1;
}
