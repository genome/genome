
package Genome::Model::Command::Services::WebApp::Loader;

use base qw( Plack::Loader::Restarter );
use Data::Dumper;

sub load {
    my $self = shift;
    my $server = shift;
    my (@args) = @_;
    $server->new(@args);
}

1;
