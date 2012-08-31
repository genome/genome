
package Bio::DNA::CoverageGap;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::CoverageGap",
    properties    => ['sequence'],
    id_properties => ['sequence'];

