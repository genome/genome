
use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 4;
use above "Genome";

use_ok('Genome::Config::Command::Grep');

my $cmd = Genome::Config::Command::Grep->create(
    query => '127786607',
    menu_item => Genome::Config::AnalysisMenu::Item->create(),
);

ok($cmd, 'created grep command');
ok($cmd->execute, 'executed grep command');

my @matching_items = $cmd->matching_items;
is(0, scalar(@matching_items), 'nothing matched search against newly minted menu item');
