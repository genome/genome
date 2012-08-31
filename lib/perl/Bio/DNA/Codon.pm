package Bio::DNA::Codon;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class
    class_name => "Bio::DNA::Codon",
    properties    => ['name'],
    id_properties => ['name'];
