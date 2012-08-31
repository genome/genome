
package Bio::DNA::Base;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::Base",
    properties    => ['name'],
    id_properties => ['name'];

