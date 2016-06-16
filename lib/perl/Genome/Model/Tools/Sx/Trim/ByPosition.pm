package Genome::Model::Tools::Sx::Trim::ByPosition;

use strict;
use warnings;

use Genome;

use Genome::Model::Tools::Sx::Split::ByNs;
use Params::Validate (qw/ :types validate_pos /);
use Regexp::Common;
use Try::Tiny;

class Genome::Model::Tools::Sx::Trim::ByPosition {
    is => 'Genome::Model::Tools::Sx::Base',
    has => {
        positions_path => {
            is => 'File',
            doc => 'Path to the file of positions to trim.',
        },    
    },
    has_transient_optional => {
        trim_positions => { is => 'HASH', },
    },
};

sub help_brief { "Trim sequences by positions" }
sub help_detail {
    <<HELP;
    The positions file is formatted with one sequence with multiple positions to trim per line. If no positions are given for a sequence, the sequence will be removed. A space follows the sequence name then the positions. Start and stop separated by a dash. Separate multiple trim positions by commas.

    Examples:

    # Trim (remove) entire SEQ1
    SEQ1

    # Trim bases 1000 to 2000 from SEQ2
    SEQ2 1000-2000

    # Trim several positions from SEQ3
    SEQ3 100-200,400-500,600-800

HELP
}

sub execute {
    my $self = shift;

    my $trim_positions = $self->load_positions( $self->positions_path );
    $self->trim_positions($trim_positions);

    $self->_init;

    my ($seqs, @seqs_to_write);
    while ( $seqs = $self->_reader->read ) {
        for my $seq ( @$seqs ) {
            $self->trim_sequence($seq);
            push @seqs_to_write, $seq if length($seq->{seq});
        }
        next if not @seqs_to_write;
        $self->_writer->write(\@seqs_to_write);
        @seqs_to_write = ();
    }

    return 1;
}

sub load_positions {
    my ($class, $positions_path) = validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});

    my $fh = Genome::Sys->open_file_for_reading($positions_path);
    my %trim_positions;
    while ( my $line = $fh->getline ) {
        chomp $line;
        my ($seq_id, $positions) = split(/\s+/, $line, 2);
        $class->fatal_message('Duplicate sequence id in trim positions! %s', $seq_id) if $trim_positions{$seq_id};

        if ( not $positions ) {
            $trim_positions{$seq_id} = [];
            next;
        }

        for my $set ( split(/,/, $positions) ) {
            my ($start, $stop) = split(/\-/, $set);
            if ( $RE{num}{int}->matches($stop) ) {
                $class->fatal_message('Invalid positions for %s! %s', $seq_id, $set) if $start == $stop or $start > $stop;
            } elsif ( $stop ne 'end' ) {
                $class->fatal_message('Invalid stop position for %s! %s', $seq_id, $stop);
            }
            push @{$trim_positions{$seq_id}}, [ $start, $stop ];
        }
    }
    $fh->close;

    return \%trim_positions;
}

sub keep_positions_for_sequence {
    my ($self, $seq) = @_;

    # Look up by id
    my $trim_positions = $self->trim_positions->{ $seq->{id} };
    if ( not $trim_positions ) {
        # Look up by orig seq id
        $trim_positions = $self->trim_positions->{ $seq->{orig_seq_id} } if $seq->{orig_seq_id};
        if ( not $trim_positions ) {
            # pcap naming
            my $pcap_seq_id = $seq->{id};
            $pcap_seq_id =~ s/scaffold/Contig/g;
            $trim_positions = $self->trim_positions->{$pcap_seq_id};
            return [ [ 0, length($seq->{seq})] ] if not $trim_positions; # keep it!
        }
    }

    return if not @$trim_positions; # sequence indicated, but no trim positions given means trim all/keep none

    my @keep_positions;
    my $current_pos = 0;
    for my $set ( @$trim_positions ) {
        push @keep_positions, [ $current_pos, ( $set->[0] - $current_pos - 1 ) ];
        $current_pos = $set->[1];
    }

    # keep the last part
    push @keep_positions, [ $current_pos, ( length($seq->{seq}) - $current_pos ) ] if length($seq->{seq}) - $current_pos > 0;
    return \@keep_positions;
}

sub trim_sequence {
    my ($self, $seq) = @_;

    my $keep_positions = $self->keep_positions_for_sequence($seq);
    if ( not $keep_positions ) { # remove
        $seq->{seq} = '';
        $seq->{qual} = '' if $seq->{qual};
        return 1;
    }

    my ($bases, $quals);
    my $has_qual = exists $seq->{qual};
    for my $set ( @$keep_positions ) {
        $bases .= substr($seq->{seq}, $set->[0], $set->[1]);
        $quals .= substr($seq->{qual}, $set->[0], $set->[1]) if $has_qual;
    }

    $seq->{seq} = $bases;
    $seq->{qual} = $quals if $has_qual;

    return 1;
}

1;

