#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Test::More tests => 12;

my @executables = Genome::Model::Tools::Cufflinks->_executable_names;

for my $app (@executables) {
    my $method = $app . '_path';

    my $v1 = Genome::Model::Tools::Cufflinks->create(use_version => '1.3.0');
    my $p1 = $v1->$method;
    is($p1, "/usr/bin/${app}1.3.0", "finds cufflinks $app path $p1 for dpkg installed cufflinks");

    my $v2 = Genome::Model::Tools::Cufflinks->create(use_version => '1.2.1');
    my $p2 = $v2->$method;
    is($p2, $ENV{GENOME_SW} . "/cufflinks/cufflinks-1.2.1.Linux_x86_64/$app", "finds $app path $p2 for legacy cufflinks installation");
}
