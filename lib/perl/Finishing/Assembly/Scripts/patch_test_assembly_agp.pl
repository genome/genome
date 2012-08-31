#!/usr/bin/env genome-perl

use strict;
use warnings;
use Say 'say';
use Data::Dumper;
use Finishing::Assembly::Factory;

my @order = (1599, 1999, 2399, 2799, 3199);

my $factory = Finishing::Assembly::Factory->connect('cmap_admin');
my $org = $factory->get_organism('test_organism');
my $chrom = $org->get_chromosome('X');
my $ass = $org->get_assembly('test_assembly');
my $scaf = $chrom->first_scaffold_for_assembly($ass)->scaffold;
my $contig = $scaf->first_contig;

say $contig->name;
$factory->schema->txn_do(
    sub{
        while (@order){
            my $scaffold_num = shift @order;
            until (! $contig->right_contig){
                $contig = $contig->right_contig;
                say $contig->name;
            }
            my $next_scaf = $ass->get_scaffold($scaffold_num);
            my $next_contig = $next_scaf->ordered_contigs->first;
            $contig->create_right_gap(dummy=>1);
            $contig->set_right_contig($next_contig);
            
        }
    }
);


=pod

=head1 NAME
ScriptTemplate - template for new perl script

=head1 SYNOPSIS

=head1 DESCRIPTION 

=cut

#$HeadURL$
#$Id$


