
package Genome::Model::Command::Services::WebApp::Resource;

use Plack::Builder;
use above 'Genome';

my $res_path = Genome::Model::Command::Services::WebApp->res_path;

builder {
    enable "Plack::Middleware::Static",
      path => qr{},
      root => $res_path;

    $app;
};
