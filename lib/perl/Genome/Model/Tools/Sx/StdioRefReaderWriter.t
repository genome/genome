#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
require IO::File;
use Test::More;

use_ok('Genome::Model::Tools::Sx::StdoutRefWriter') or die;
is(Genome::Model::Tools::Sx::StdoutRefWriter->type, 'ref', 'type is ref');
my $writer = Genome::Model::Tools::Sx::StdoutRefWriter->create();
ok($writer, 'Created writer');
can_ok($writer, 'write');

use_ok('Genome::Model::Tools::Sx::StdinRefReader') or die;
is(Genome::Model::Tools::Sx::StdinRefReader->type, 'ref', 'type is ref');
my $reader = Genome::Model::Tools::Sx::StdinRefReader->create();
ok($reader, 'Created reader');
can_ok($reader, 'read');

my $fh = IO::File->new(qq{ $^X -Mstrict -M'above "Genome"' -MGenome::Model::Tools::Sx::StdoutRefWriter  -e 'Genome::Model::Tools::Sx::StdoutRefWriter->write(UR::Value->get(100));' |  $^X -Mstrict -M'above "Genome"' -MGenome::Model::Tools::Sx::StdinRefReader -e 'my \$ref = Genome::Model::Tools::Sx::StdinRefReader->read or die; print \$ref->id."\n";' | });
ok($fh, 'Created pipe') or die;
my $value = $fh->getline;
$fh->close;
chomp $value;
is($value, 100, 'Write object to STDIN, read object from STDOUT, got correct value');

done_testing();
