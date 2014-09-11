package Genome::File::Breakdancer::ToBedPe;

use strict;
use warnings;

use Carp qw(confess);
use Genome::File::BedPe::Entry;

sub new {
    my ($class, $header, $slop, $name_generator) = @_;

    my $self = {
        _header => $header,
        _slop => $slop,
        _name_generator => $name_generator,
        _generated => {},
    };

    bless $self, $class;
    return $self;
}

sub add_generated {
    my ($self, $name, $code) = @_;
    if (exists $self->{_generated}{$name}) {
        confess "A generated field with name $name already exists!"
    }

    $self->{_generated}{$name} = $code;
}

sub convert {
    my ($self, $entry) = @_;

    my $slop = $self->{_slop};
    my $start1 = $entry->{pos1} > $slop ? $entry->{pos1} - $slop : 0;
    my $start2 = $entry->{pos2} > $slop ? $entry->{pos2} - $slop : 0;
    my %fields = (
        chrom1 => $entry->{chr1},
        start1 => $start1,
        end1 => $entry->{pos1} + $self->{_slop},
        chrom2 => $entry->{chr2},
        start2 => $start2,
        end2 => $entry->{pos2} + $self->{_slop},
        score => $entry->{score}
    );

    my $name = $self->{_name_generator}->($entry);
    if ($name) {
        $fields{name} = $name;
    }

    my $gen = $self->{_generated};
    for my $name (keys %$gen) {
        my $value = $gen->{$name}->($entry);
        $fields{$name} = $value;
    }

    return Genome::File::BedPe::Entry->new_from_fields($self->{_header}, %fields);
}

1;
