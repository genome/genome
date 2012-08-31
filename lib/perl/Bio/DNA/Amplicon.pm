
package Bio::DNA::Amplicon;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class
    class_name => "Bio::DNA::Amplicon",
    properties    => ['primer1','primer2'],    
    id_properties => ['primer1','primer2'];

