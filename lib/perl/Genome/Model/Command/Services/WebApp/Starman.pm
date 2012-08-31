
package Genome::Model::Command::Services::WebApp::Starman;

use strict;
use warnings;

use Genome::Model::Command::Services::WebApp::Starman::Server;

sub new {
    my $class = shift;
    bless { @_ }, $class;
}

sub run {
    my($self, $app) = @_;

    if ($ENV{SERVER_STARTER_PORT}) {
        require Net::Server::SS::PreFork;
        @Starman::Server::ISA = qw(Net::Server::SS::PreFork); # Yikes.
    }

    Genome::Model::Command::Services::WebApp::Starman::Server->new->run($app, {%$self});
}

1;
