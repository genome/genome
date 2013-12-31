
#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More;

my $class = 'Genome::Config::AnalysisProject::Command::ShowConfig';
use_ok($class);

done_testing();

