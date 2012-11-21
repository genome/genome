#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Site::TGI::Synchronize::ReconcileMiscUpdate') or die;

# TODO Create entities
# TODO Create misc updates

my $reconcile = Genome::Site::TGI::Synchronize::ReconcileMiscUpdate->create();
#ok($reconcile, 'create reconcile command');
#ok($reconcile->execute, 'create reconcile command');

done_testing();
exit;

