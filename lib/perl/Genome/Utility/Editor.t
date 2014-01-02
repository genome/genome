#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above 'Genome';

my $editor = $ENV{'EDITOR'};
$ENV{'EDITOR'} = 'touch';

my $class = 'Genome::Utility::Editor';
use_ok($class);

is(Genome::Utility::Editor::from_empty_file(), 0,
    'expect 0 return code for empty file');

is(Genome::Utility::Editor::from_empty_file(allow_empty => 1), '',
    'expect empty string when allow_empty is on for empty file');

my $existing_file = Genome::Sys->create_temp_file_path();
my $existing_contents = 'This string exists!';

Genome::Sys->write_file($existing_file, $existing_contents);

is(Genome::Utility::Editor::from_existing_contents($existing_contents), $existing_contents,
    'expect existing contents to come back identical');

is(Genome::Utility::Editor::from_existing_file($existing_file), $existing_contents,
    'expect the the sub to return the contents');

$ENV{'EDITOR'} = $editor;

done_testing();
