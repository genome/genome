
package Bio::DNA::CombinedCDSRegion;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::CombinedCDSRegion",
    properties    => ['name'],
    id_properties => ['name'];

