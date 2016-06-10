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
            push @{$trim_positions{$seq_id}}, [qw/ 1 end /];
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

1;

