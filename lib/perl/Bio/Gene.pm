
package Bio::Gene;
use base 'Bio::Dynamic';

use strict;
use warnings;

create App::Object::Class
    class_name => "Bio::Gene",
    properties    => ['locus_link'],
    id_properties => ['locus_link'];

