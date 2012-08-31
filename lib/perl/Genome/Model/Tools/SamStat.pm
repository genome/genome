package Genome::Model::Tools::SamStat;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SamStat {
    is  => 'Command::Tree',
    doc => 'Displaying sequence statistics for next generation sequencing',
};

sub help_brief {
    "Run tools and handle output from samstat: http://samstat.sourceforge.net/",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt sam-stat ...    
EOS
}
sub help_detail {
    "For each input file SAMStat will create a single html page named after the input file name plus a dot html suffix. http://samstat.sourceforge.net/";
}


1;
