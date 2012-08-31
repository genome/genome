use warnings;
use strict;

use File::Basename;
use lib File::Basename::dirname(__FILE__)."/../../../lib";
use lib File::Basename::dirname(__FILE__)."/../..";
use UR;
use Test::More tests => 4;

UR::Object::Type->define(
    class_name => 'Test::Value',
    is => 'UR::Value',
    id_by => [
        string => { is => 'Text' }
    ]
);

my $z = Test::Value->get('xyz');
ok($z,"get('xyz') returned on first call");

my $y = Test::Value->get('xyz');
ok($y,"get('xyz') returned on second call");

my $x = Test::Value->get(string => 'abc');
ok($x,"get(string => 'abc') returned on first call");

my $w = Test::Value->get(string => 'abc');
ok($w,"get(string => 'abc') returned on second call");

 
