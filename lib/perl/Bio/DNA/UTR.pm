
package Bio::DNA::UTR;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class    
    class_name => "Bio::DNA::UTR",
    properties    => ['name'],
    id_properties => ['name'];

