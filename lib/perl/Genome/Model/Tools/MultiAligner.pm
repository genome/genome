package Genome::Model::Tools::MultiAligner;

use strict;
use warnings;

use Genome;
use File::Basename;

my $DEFAULT = '0';

class Genome::Model::Tools::MultiAligner {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "No version needed, ignore this" },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run multiple aligners",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools multialigner ...    
EOS
}

sub help_detail {                           
    return <<EOS 
FIXME
EOS
}

sub path_for_multialigner_version {
    return '';
}

sub default_multialigner_version {
    return $DEFAULT;
}
        
sub default_version { 
    return $DEFAULT; 
}

1;

