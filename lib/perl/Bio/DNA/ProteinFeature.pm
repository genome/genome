
package Bio::DNA::ProteinFeature;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::ProteinFeature",
    properties    => ['name'],
    id_properties => ['name'];

