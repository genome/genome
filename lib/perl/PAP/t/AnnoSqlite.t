use strict;
use warnings;

use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use Cwd;
use File::Temp;
use Test::More tests => 2;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::AnnoSqlite');
}


# need to try to run annosqlite...
