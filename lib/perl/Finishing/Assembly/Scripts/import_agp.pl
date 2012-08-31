#!/usr/bin/env genome-perl
use strict;
use warnings;
use Finishing::Assembly::AGP::Importer;
use Finishing::Assembly::Factory;
use Finishing::Assembly::AGP::Reader;
my $factory = Finishing::Assembly::Factory->connect('cmap_admin');  
my $org_name = 'test_organism'; 
my $assembly_name = 'test_assembly';
my $file = '/gscuser/adukes/bin/testers/AssemblyImprovement/Source/chimp_050412.pcap.contigs399.agp';

my $reader = Finishing::Assembly::AGP::Reader->new(io => $file);
my $organism = $factory->get_organism('test_organism');
my $assembly = $organism->get_assembly($assembly_name);

my $importer = Finishing::Assembly::AGP::Importer->new(
    reader => $reader, 
    assembly => $assembly,
    txn_do => sub { $assembly->schema->txn_do(@_)},
    ordered => 1, 
);

$importer->execute;

