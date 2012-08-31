package Genome::Model::Tools::SeeFive::Rule;

use strict;
use warnings;
use Genome;            

class Genome::Model::Tools::SeeFive::Rule {
    is => 'Command',
    has => [ 
        trial           => { is => 'Genome::Model::Tools::SeeFive::Rule', id_by => 'trial_id' },
        n               => {},
        possible_items  => {},
        wrong_items     => {},
        lift            => {},
        lines           => {},
    ],
};

sub perl_src {
    my $self = shift;
    my $lines = $self->lines;
    chomp @$lines;
    unless (@$lines) {
        return "# no lines $self->{n}?";
    }
    my $src;
    my %metric_names;
    for my $line (@$lines) {
        if ($line =~ /\s+\-\>/) {
            my ($class,$prob) = ($line =~ /\s+\-\>\s*class\s+(\S+)\s+\[([\d+\.]+)\]/);
            $src .= "\n) ? (\"$class\", $prob) : ()\n";
            return $src;
        }
        my ($metric,$op,$value) = ($line =~ /\s+([\w\-]+)\s+(\S+)\s+(\S+)/);
        $metric =~ s/\s/SPACE/g;
        $metric =~ s/\+/PLUS/g;
        $metric =~ s/\-/MINUS/g;
        unless ($metric) {
            die "failed to parse line: $line\n";
        }
        $metric_names{$metric}++;
        $src .= ($src ? "\n    and (" : "(\n    (");
        $src .= '$' . $metric . " $op $value)";
    }
    die "No final value for rule?";
}

1;

