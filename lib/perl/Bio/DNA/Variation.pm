
package Bio::DNA::Variation;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::Variation",
    properties    => ['name'],
    id_properties => ['name'];

