
package Genome::Model::Command::Services::WebApp::404Handler;

use strict;
use warnings;

sub {
    my ( $content ) = @_;

    my $string = join( "\n", @$content );
    my $doc = <<"    HTML";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
    <!--template: status/root.xsl:match "/"-->
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <title>Error: $string</title>
    <link rel="shortcut icon" href="/res/img/gc_favicon.png" type="image/png" />
    <link rel="stylesheet" href="/res/css/blueprint/screen.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/blueprint/print.css" type="text/css" media="print" />
    <link rel="stylesheet" href="/res/css/master.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/buttons.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/icons.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/forms.css" type="text/css" media="screen, projection" />
    <link type="text/css" href="/res/css/jquery-ui.css" rel="stylesheet" />
    <link href="/res/css/jquery-ui-overrides.css" type="text/css" rel="stylesheet" media="screen, projection" />
  </head>
  <body>
    <div class="page">
      <div class="header rounded-bottom gradient-grey shadow">
        <div class="container">
          <div class="title span-24 last app_error_32">
            <h1>Error Encountered</h1>
          </div>
        </div>
      </div>

      <div class="content rounded shadow" style="background-color: #FAA">
        <div class="container">
        <div class="span-24 last">
          <div class="rounded" style="background: #FFF; margin-bottom: 10px;">
            <div class="padding10">
              <p><strong>Error:</strong> $string<br>
                Please email apipe\@genome.wustl.edu</p>
            </div>
          </div>
        </div>
        </div>
      </div>
    </div>
  </body>
</html>
    HTML

    [ 404, [ 'Content-type', 'text/html' ], [$doc] ];
};
