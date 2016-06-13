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
            doc => 'Path to he file of positions to trim.',
        },    
    },
    has_transient_optional => {
        trim_positions => { is => 'HASH', },
    },
};

sub help_brief { "Trim sequences by positions" }
sub help_detail { help_brief() }

sub load_positions {
    my ($class, $positions_path) = validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});

    my $fh = Genome::Sys->open_file_for_reading($positions_path);
    my %trim_positions;
    while ( my $line = $fh->getline ) {
        chomp $line;
        my ($seq_id, $positions) = split(/\s+/, $line, 2);
        $class->fatal_message('Duplicate sequence id in trim positions! %s', $seq_id) if $trim_positions{$seq_id};
        
        if ( not $positions ) {
            $trim_positions{$seq_id} = 'ALL';
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
            return 'ALL' if not $trim_positions; # keep it!
        }
    }

    return if not ref $trim_positions; # trim all

    my @keep_positions;
    my $current_pos = 0;
    for my $set ( @$trim_positions ) {
        push @keep_positions, [ $current_pos, ( $set->[0] - $current_pos - 1 ) ];
        $current_pos = $set->[1];
    }

    # keep the last part
    push @keep_positions, [ $current_pos, ( length($seq->{seq}) - $current_pos ) ];

    return \@keep_positions;
}

1;

