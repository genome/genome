
package Genome::Model::Command::Services::WebApp::Command::Run;

class Genome::Model::Command::Services::WebApp::Command::Run {
    is => 'Genome::Command::Base',
    has_optional => [
        fixed_port => {
            is => Boolean,
            default => 0,
            doc => "force a fixed port, useful in daemon mode",
        },
    ]
};

sub execute {
  my $self = shift;
  eval {
     $app = Genome::Model::Command::Services::WebApp->create( fixed_port => $self->fixed_port );
     $app->execute();
  };
  if ($@) {
    die "Failed during execute(): $@";
  }
}

1;
