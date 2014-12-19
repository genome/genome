package Genome::File::BamReadcount::IndexedReader;

use strict;
use warnings;
use Genome::File::BamReadcount::Reader;

sub new {
    my ($class, @args) = @_;
    my $reader = Genome::File::BamReadcount::Reader->new(@args);

    my $self = {
        _reader => $reader,
        _prev_entry => undef,
        _cur_entry => undef,
        _prev_pos => -1,
        _cur_pos => -1,
        _cur_chrom => 1,
        _cur_chroms_seen => Set::Scalar->new,
    };

    bless $self, $class;
    return $self;
}

sub get_entry {
    my ($self, $chrom, $pos) = @_;

    if ($chrom ne $self->{_cur_chrom} and $self->{_cur_chroms_seen}->contains($chrom)) {
        die sprintf("Chromosome %s has already been passed by", $chrom);
    }

    while ($pos > $self->{_cur_pos} or $chrom ne $self->{_cur_chrom}) {
        $self->read_line;
        return unless defined $self->{_cur_entry};
    }

    if ($pos == $self->{_prev_pos}) {
        return $self->{_prev_entry};
    }
    elsif ($pos == $self->{_cur_pos}) {
        return $self->{_cur_entry};
    }
    elsif ($pos < $self->{_prev_pos}) {
        die sprintf("Position %s has already been passed by", $pos);
    }
    else {
        return;
    }
}

sub read_line {
    my $self = shift;
    $self->{_prev_entry} = $self->{_cur_entry};
    $self->{_prev_pos} = $self->{_cur_pos};

    $self->{_cur_entry} = $self->{_reader}->next;
    if (defined $self->{_cur_entry}) {
        if ($self->{_cur_pos} == $self->{_cur_entry}->{_position} and
            $self->{_cur_chrom} eq $self->{_cur_entry}->{_chromosome} ) {
            die sprintf("Duplicate entries detected at chromosome (%s) position (%s).
                Bam readcount file must be deduplicated",
                $self->{_cur_chrom},
                $self->{_cur_pos});
        }
        $self->{_cur_chrom} = $self->{_cur_entry}->{_chromosome};
        $self->{_cur_chroms_seen}->insert($self->{_cur_chrom});
        $self->{_cur_pos} = $self->{_cur_entry}->{_position};
    }
    return;
}

1;

