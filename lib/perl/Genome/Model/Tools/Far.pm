package Genome::Model::Tools::Far;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Far {
    is  => 'Command::Tree',
    doc => 'Flexible Adapter Remover(FAR)',
};

sub help_brief {
    "To trim adaptor sequences",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt far ...    
EOS
}

sub help_detail {
    return <<EOS 
More information about the FAR tool can be found at http://sourceforge.net/apps/mediawiki/theflexibleadap/index.php?title=Main_Page.
EOS
}

1;