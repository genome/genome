#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;
use Finishing::Assembly::Factory;
use Finishing::Assembly::Ace::Exporter;

my $factory = Finishing::Assembly::Factory->connect('cmap_admin');
my $file = '/gscuser/adukes/svn/test_modules/Finishing/Assembly/test/TestAssembly/assembly.ace';
unlink $file;


my $organism = $factory->get_organism('test_organism');
my $assembly = $organism->get_assembly('test_assembly');
my $chromosome = $organism->get_chromosome('X');
my $scaffold = $chromosome->first_scaffold_for_assembly($assembly)->scaffold;

my $contig = $scaffold->first_contig;

my $exporter = Finishing::Assembly::Ace::Exporter->new(file => $file);

$exporter->export_contig(contig => $contig);
while ($contig = $contig->right_contig and defined $contig){
    $exporter->export_contig(contig => $contig);
}

$exporter->close;


=pod

=head1 NAME
ScriptTemplate - template for new perl script

=head1 SYNOPSIS

=head1 DESCRIPTION 

=cut

#$HeadURL$
#$Id$


