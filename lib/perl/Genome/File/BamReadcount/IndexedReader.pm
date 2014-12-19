package Genome::File::BamReadcount::IndexedReader;

use strict;
use warnings;
use Genome::File::BamReadcount::Reader;

sub new {
    my ($class, @args) = @_;
    my $reader = Genome::File::BamReadcount::Reader->new(@args);

    my $self = {
        _reader => $reader,
        _previous_entry => undef,
        _current_entry => undef,
        _previous_position => -1,
        _current_position => -1,
        _current_chromosome => 1,
        _chromosomes_seen => Set::Scalar->new,
    };

    bless $self, $class;
    return $self;
}

sub get_entry {
    my ($self, $chrom, $pos) = @_;

    if ($chrom ne $self->{_current_chromosome} and $self->{_chromosomes_seen}->contains($chrom)) {
        die sprintf("Chromosome %s has already been passed by", $chrom);
    }

    while ($pos > $self->{_current_position} or $chrom ne $self->{_current_chromosome}) {
        $self->read_line;
        return unless defined $self->{_current_entry};
    }

    if ($pos == $self->{_previous_position}) {
        return $self->{_previous_entry};
    }
    elsif ($pos == $self->{_current_position}) {
        return $self->{_current_entry};
    }
    elsif ($pos < $self->{_previous_position}) {
        die sprintf("Position %s has already been passed by", $pos);
    }
    else {
        return;
    }
}

sub read_line {
    my $self = shift;
    $self->{_previous_entry} = $self->{_current_entry};
    $self->{_previous_position} = $self->{_current_position};

    $self->{_current_entry} = $self->{_reader}->next;
    if (defined $self->{_current_entry}) {
        if ($self->{_current_position} == $self->{_current_entry}->{_position} and
            $self->{_current_chromosome} eq $self->{_current_entry}->{_chromosome} ) {
            die sprintf("Duplicate entries detected at chromosme (%s) position (%s).
                Bam readcount file must be deduplicated",
                $self->{_current_chromosome},
                $self->{_current_position});
        }
        $self->{_current_chromosome} = $self->{_current_entry}->{_chromosome};
        $self->{_chromosomes_seen}->insert($self->{_current_chromosome});
        $self->{_current_position} = $self->{_current_entry}->{_position};
    }
    return;
}

1;

