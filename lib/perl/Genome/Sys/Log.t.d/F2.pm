use strict;
use warnings;
package F2;
sub f2 {
    my ($type,$msg) = @_;
    my $method = $type . "_message";
    Genome->$method($msg);
}
1;

