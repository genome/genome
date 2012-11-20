package Genome::Model::Tools::Annovar::Base;

use strict;
use warnings;
use Genome;

my ($ANNOVAR_DIR) = Cwd::abs_path(__FILE__) =~ /(.*)\//;
my $ANNOVAR_SCRIPT_PATH = $ANNOVAR_DIR . "/Annovar.d/";

class Genome::Model::Tools::Annovar::Base {
    is => "Command::V2",
    has => [
    ],
};

sub script_path {
    return $ANNOVAR_SCRIPT_PATH;
}

1;

