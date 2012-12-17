#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;
#use Test::More tests => 6;
use Test::More tests => 1;


BEGIN {
    use_ok('Genome::Model::Tools::Hgmi::Hap');
}

#my $testrunconfig = "/gsc/var/cache/data/Genome-Model-Tools-Hgmi/hap-test0.yaml";
#my $h = Genome::Model::Tools::Hgmi::Hap->create( gen_example => 1,
#                                                 config => '/tmp/testconfig.yaml',
#                                                 dev => 1);
#
#ok($h); # isa_ok???
#ok($h->execute());
#ok(-f '/tmp/testconfig.yaml');
#remove \1, qw{ /tmp/testconfig.yaml };
#
#my $h2 = Genome::Model::Tools::Hgmi::Hap->create( 
#              config => "/tmp/nonexistconfig.yaml",
#              dev => 1 );
#
#ok($h2);
#is($h2->execute(), 0,'complain about non-existent config file');
#
#my $h3 = Genome::Model::Tools::Hgmi::Hap->create( config => $testrunconfig,
#                                                  dev => 1 );
#
#ok($h3->execute());
##

1;
