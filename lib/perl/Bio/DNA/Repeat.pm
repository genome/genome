
package Bio::DNA::Repeat;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::Repeat",
    properties    => ['name'],
    id_properties => ['name'];

