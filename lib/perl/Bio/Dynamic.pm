
package Bio::Dynamic;
use base 'App::Object';

use strict;
use warnings;

define App::Object::Class
    class_name => "Bio::Dynamic",
    properties    => ['obj_id'],
    id_properties => ['obj_id'];

sub load {
    my $class = shift;
    my $params = $class->preprocess_params(@_);
    
    # See if the requested object is loaded.
    my @loaded = $class->is_loaded($params);
    $class->context_return(@loaded) if @loaded;
    
    # If there is no ID, we can't autogenerate.
    unless ($params->{id}) {
        return ();
    }
    
    # Auto generate the object on the fly.
    my $obj = $class->create_object($params);
    $obj->signal_change("load");
    return $obj;
}

1;