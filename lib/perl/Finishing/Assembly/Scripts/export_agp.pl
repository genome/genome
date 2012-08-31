#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;
use Finishing::Assembly::Factory;
use Finishing::Assembly::AGP::Exporter;

my $factory = Finishing::Assembly::Factory->connect('cmap_admin');
my $assembly = $factory->get_organism('test_organism')->get_assembly('test_assembly');

my $io = IO::Handle->new;
$io->fdopen(fileno(STDOUT), 'w');
die unless $io;


my $exporter = Finishing::Assembly::AGP::Exporter->new(
    ordered=>1,
    assembly=>$assembly,
    chromosome_name=>'X',
    writer=>Finishing::Assembly::AGP::Writer->new(io => $io),
);

$exporter->execute;

=pod

=head1 NAME
ScriptTemplate - template for new perl script

=head1 SYNOPSIS

=head1 DESCRIPTION 

=cut

#$HeadURL$
#$Id$


