
package Bio::DNA::CombinedExonRegion;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::CombinedExonRegion",
    properties    => ['name'],
    id_properties => ['name'];

