package Genome::File::BamReadcount::BufferedReader;

use strict;
use warnings;
use Genome::File::BamReadcount::Reader;
use Genome::File::Location;
use Data::Dump qw(pp);

sub new {
    my ($class, @args) = @_;
    my $reader = Genome::File::BamReadcount::Reader->new(@args);

    my $self = {
        _reader => $reader,
        _prev_entry => undef,
        _cur_entry => undef,
        _prev_loc => Genome::File::Location->new(-1, -1),
        _cur_loc => Genome::File::Location->new(-1, -1),
    };

    bless $self, $class;
    return $self;
}

sub get_entry {
    my ($self, $chrom, $pos) = @_;

    my $target_loc = Genome::File::Location->new($chrom, $pos);
    while ($target_loc->is_greater_than($self->{_cur_loc})) {
        $self->read_line;
        return unless defined $self->{_cur_entry};
    }

    if ($target_loc->is_equal_to($self->{_prev_loc})) {
        return $self->{_prev_entry};
    }
    elsif ($target_loc->is_equal_to($self->{_cur_loc})) {
        return $self->{_cur_entry};
    }
    elsif ($target_loc->is_less_than($self->{_prev_loc})) {
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

    $params{file_name} = $params{_reader}->{name};
    delete $params{_reader};

    $params{_prev_loc} = $params{_prev_loc}->to_string();
    $params{_cur_loc} = $params{_cur_loc}->to_string();
    return pp(\%params);
}

sub read_line {
    my $self = shift;
    $self->{_prev_entry} = $self->{_cur_entry};
    $self->{_prev_loc}->set_equal_to($self->{_cur_loc});

    my $entry = $self->{_reader}->next;
    $self->{_cur_entry} = $entry;
    if (defined $entry) {
        my $entry_loc = Genome::File::Location->new(
            $entry->{_chromosome},
            $entry->{_position},
        );
        if ($self->{_cur_loc}->is_equal_to($entry_loc)) {
            die sprintf("Duplicate entries detected at %s.  " .
                "Bam readcount file (%s) must be deduplicated",
                $self->{_cur_loc},
                $self->{_reader}->{name}
            );
        } else {
            $self->{_cur_loc}->set_equal_to($entry_loc);
        }
    }
    return;
}

1;

