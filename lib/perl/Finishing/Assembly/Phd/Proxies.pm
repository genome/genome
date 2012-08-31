package Finishing::Assembly::Phd::AssembledReadProxy;

use base 'Finishing::Assembly::Proxy';

use Data::Dumper;

sub methods_for_source_method : RESTRICTED
{
    return 
    (
        name => undef,
        base_string => undef,
        qualities => undef,
        position => undef,
        complemented => undef,
        tags => 'read_tag',
        'time' => undef,
        chromat_file => undef,
        phd_file => undef,
        chem => undef,
        dye => undef,
        'length' => undef,
    );
}

1;

#$HeadURL$
#$Id$
