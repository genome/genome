#!/usr/bin/env genome-perl
use strict;
use warnings;

use Test::More;
use above "Genome";
use Genome::Utility::Test;

my $class = 'Genome::Config::Parser::YAML';

use_ok($class);

my $data_dir = Genome::Utility::Test->data_dir($class, 1);
ok(-d $data_dir, "data_dir exists: $data_dir") or die;

sub parse_in_eval {
    my ($filename, $message) = @_;
    eval {
        Genome::Config::Parser::YAML->parse($data_dir . '/' . $filename);
    };
    ok($@, $message);
}

sub parse_valid {
    my ($filename, $message) = @_;
    ok(Genome::Config::Parser::YAML->parse($data_dir . '/' . $filename), $message);
}

parse_in_eval('other','Parser will die when given a file without a valid YAML extension');
parse_in_eval('empty.yml','Parser will die when given an empty file');

parse_valid('valid.yml', 'Parser accepts upper and lower case extensions');
parse_valid('valid.YML', 'Parser accepts upper and lower case extensions');
parse_valid('valid.yaml', 'Parser accepts both extensions');

done_testing();
