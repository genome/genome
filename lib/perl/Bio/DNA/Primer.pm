
package Bio::DNA::Primer;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class
    class_name => "Bio::DNA::Primer",
    properties    => ['sequence'],
    id_properties => ['sequence'];
