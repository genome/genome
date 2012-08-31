package Genome::Model::Tools::RefCov::ROI::Bed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::ROI::Bed {
    is => ['Genome::Model::Tools::RefCov::ROI::FileI'],
};

sub _load_chromosomes {
    my $self = shift;
    my $cmd = 'cut -f 1 '. $self->file .' | sort -u';
    my $output = `$cmd`;
    my @chromosomes = split("\n",$output);
    $self->_chromosomes(\@chromosomes);
}

sub _parse_line {
    my $self = shift;
    my $line = shift;
    chomp($line);
    my ($chr,$start,$end,$name,$score,$strand) = split("\t",$line);
    unless (defined($chr) && defined($start) && defined($end)) {
        return;
    }
    #BED format uses zero-based start coordinate, convert to 1-based
    $start += 1;
    my $wingspan = $self->wingspan;
    if ($wingspan) {
        $start -= $wingspan;
        $end += $wingspan;
    }
    if ( !defined($name) || $name eq '') {
        $name = $chr .':'. $start .'-'. $end;
    }
    my %region = (
        name => $name,
        chrom => $chr,
        start => $start,
        end => $end,
    );
    if (defined($strand)) {
        $region{strand} = $strand;
    }
    return \%region;
}

1;
