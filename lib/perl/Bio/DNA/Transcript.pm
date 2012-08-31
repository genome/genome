
package Bio::DNA::Transcript;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::Transcript",
    properties    => ['name'],
    id_properties => ['name'];

