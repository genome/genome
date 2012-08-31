package Genome::InstrumentData::Command::Align::Clc;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::Command::Align::Clc {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name => { value => 'clc' },
    ],
    has_param => [
        version      => { default_value => '3.11.50018' },
    ],
    doc => "align instrument data using ClcBio assembler/aligner"
};

sub help_synopsis {
    return <<EOS
This runs the Clc aligner in either paired end mode or fragment mode (depending on how many inputs are fed in).  This module uses hardcoded insert size values of 180 - 250bp.   Original output is produced in cas format, and then converted into (headerless) sam format.

EOS
}

sub help_detail {
return <<EOS
Launch the Clc aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}

1;

