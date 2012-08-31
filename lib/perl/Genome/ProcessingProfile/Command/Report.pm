package Genome::ProcessingProfile::Command::Report;

#REVIEW fdu 11/20/2009
#1. Remove 'use Data::Dumper';
#2. Add help info
#3. This module can be dropped (see notes on G::P::C::R::Summary)

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::ProcessingProfile::Command::Report {
    is => 'Command',
};

1;

