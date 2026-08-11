package Genome::Sys::SLURM::JobIterator;

use strict;
use warnings;

package Job::Iterator;

sub new {
    my $class = shift;
    my $opts = shift || '';
    open my $f, "squeue -O '" . join(":\t,", qw(JOBID PARTITION ACCOUNT USERNAME STATE)) . ":' " . $opts . ' |';
    my $headers = readline($f);
    chomp $headers;
    return bless[ $f, [split "\t", $headers]], $class;
}

sub next {
    my $self = shift;

    my $line = readline($self->[0]);
    return unless defined $line;

    chomp $line;
    return Job->new($self->[1], $line);
}

package Job;

sub new {
    my $class = shift;
    my $headers = shift;
    my $line = shift;

    my @fields = split("\t", $line);

    my %data;
    for my $i (0..(scalar(@fields)-1)) {
        $data{ $headers->[$i] } = $fields[$i];
    }

    #for compatibility
    $data{Job} = $data{JOBID};
    $data{Status} = $data{STATE};

    return bless \%data, $class;
}

sub started_on {
    die('started_on unimplemented for slurm iterator');
}

sub events {
    die('events unimplemented for slurm iterator');
}

1;
