
package Bio::DNA::Clone;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::Clone",
    properties    => ['name'],
    id_properties => ['name'];

