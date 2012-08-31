
package Bio::DNA::Exon;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::Exon",
    properties    => ['name'],
    id_properties => ['name'];

