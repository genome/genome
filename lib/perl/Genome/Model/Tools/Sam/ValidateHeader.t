#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

BEGIN {
    if (`uname -a` =~ /x86_64/){
        plan tests => 5;
    }
    else{
        plan skip_all => 'Must run on a 64 bit machine';
    }

    use_ok('Genome::Model::Tools::Sam::ValidateHeader');
}

my $input = Genome::Config::get('test_inputs') . "/Genome-InstrumentData-Alignment/new.bam";
my $input_bad = Genome::Config::get('test_inputs') . "/Genome-InstrumentData-Alignment/old.bam";

my $hv = Genome::Model::Tools::Sam::ValidateHeader->create(input_file=>$input);
ok($hv, "created command");

my $result = $hv->execute;
ok($result, "executed, header ok");

my $hv2 = Genome::Model::Tools::Sam::ValidateHeader->create(input_file=>$input_bad);
ok($hv2, "created command");

my $result2 = $hv2->execute;
ok(!defined($result2), "executed, successfully reported bad header");
