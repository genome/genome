
package Genome::Model::Command::Services::WebApp::Runner;

use base qw( Plack::Runner );
use Genome::Model::Command::Services::WebApp::Starman;

sub load_server {
    my($self, $loader) = @_;
    $self->{server}->new(@{$self->{options}});
}

1;
