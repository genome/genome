package Genome::File::BamReadcount::BufferedReader;

use strict;
use warnings;
use Genome::File::BamReadcount::Reader;
use Set::Scalar;
use Sort::strverscmp;

sub new {
    my ($class, @args) = @_;
    my $reader = Genome::File::BamReadcount::Reader->new(@args);

    my $self = {
        _reader => $reader,
        _prev_entry => undef,
        _cur_entry => undef,
        _prev_pos => -1,
        _prev_chrom => -1,
        _cur_pos => -1,
        _cur_chrom => 1,
    };

    bless $self, $class;
    return $self;
}

sub get_entry {
    my ($self, $chrom, $pos) = @_;

    while (
        ($chrom eq $self->{_cur_chrom} and $pos > $self->{_cur_pos}) or
        strverscmp($chrom, $self->{_cur_chrom}) == 1
    ) {
        $self->read_line;
        return unless defined $self->{_cur_entry};
    }
    if ($pos == $self->{_prev_pos} and $chrom eq $self->{_prev_chrom}) {
        return $self->{_prev_entry};
    }
    elsif ($pos == $self->{_cur_pos} and $chrom eq $self->{_cur_chrom}) {
        return $self->{_cur_entry};
    }
    elsif (
        ($chrom eq $self->{_prev_chrom} and $pos < $self->{_prev_pos}) or
        strverscmp($chrom, $self->{_prev_chrom}) == -1
    ) {
        die sprintf("Position %s on chromosome %s has already been passed by.  Current state: %s",
            $pos,
            $chrom,
            $self->to_string);
    }
    else {
        return;
    }
}

sub to_string {
    my $self = shift;
    my %params = %$self;
    delete $params{_cur_entry};
    delete $params{_prev_entry};
    delete $params{_reader};
    $params{file_name} = $self->{_reader}->{name};
    Data::Dumper::Dumper(\%params);
}

sub read_line {
    my $self = shift;
    $self->{_prev_entry} = $self->{_cur_entry};
    $self->{_prev_pos} = $self->{_cur_pos};
    $self->{_prev_chrom} = $self->{_cur_chrom};

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
        $self->{_cur_pos} = $self->{_cur_entry}->{_position};
    }
    return;
}

1;

