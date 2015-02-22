package Genome::InstrumentData::Command::Remove;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Remove {
    is => 'Command::V2',
    has_input => [
        instrument_data => { 
            is => 'Genome::InstrumentData::Imported',
            shell_args_position => 1,
            is_many => 1,
            doc => 'instrument data to remove, specified by id or expression'
        },
    ],
    doc => 'remove imported instrument data',
};

sub execute {
    my $self = shift;
    my @i = $self->instrument_data();
    $self->status_message("Removing " . scalar(@i) . " instrument data entries...");
    sleep 5;
    for my $i (@i) {
        $self->status_message("deleting " . $i->__display_name__ . "...");
        $i->delete;
        print "$i\n";
    }
    return 1;
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/InstrumentData/Command/Remove.pm $
#$Id: Remove.pm 53285 2009-11-20 21:28:55Z fdu $
