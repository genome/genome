#!/usr/bin/env genome-perl

use strict;
use warnings;


use above "PAP";
use above 'Workflow';
use Test::More qw(no_plan);
use Data::Dumper;
use PAP;

ok(1);
# this needs to be updated.
# these objects have changed.
#my $i = Workflow::Store::Db::Operation::Instance->get(88);
#ok($i->operation->set_all_executor(Workflow::Executor::SerialDeferred->create()));

#$i->status('crashed');
#$i->is_done(0);

#my $mainpeer = Workflow::Store::Db::Operation::Instance->get(92);
#for my $p ($mainpeer, $mainpeer->peers) {
#    $p->is_done(0);
#    $p->status('new');
#    $p->output_connector->status('new');
#    $p->output_connector->is_done(0);
#}

#my $dbupload = Workflow::Store::Db::Operation::Instance->get(101);
#$dbupload->status('crashed');
#$dbupload->is_done(0);

#for my $id (90,99,89,100,91) {
#    my $obj = Workflow::Store::Db::Operation::Instance->get($id);

#    $obj->status('new');
#    $obj->is_done(0);
#}

#my $outcat = Workflow::Store::Db::Operation::Instance->get(90);
#my $mainpeer = Workflow::Store::Db::Operation::Instance->get(92);
#foreach my $p ($mainpeer, $mainpeer->peers) {
#    $outcat->input()->{'blastp psortb feature'}->[$p->parallel_index] =
#    Workflow::Link::Instance->create(
#        operation_instance => $p,
#        property => 'bio seq features'
#    );
#}

#ok($i->resume());
#ok($i->operation->wait);

#done_testing();
1;
