#!/usr/bin/env genome-perl

use strict;
use warnings;
use Finishing::Assembly::Factory;

my $factory = Finishing::Assembly::Factory->connect('cmap_admin');

my $organism = $factory->get_organism("pan_troglodytes");

my $ai = $organism->assemblies;

my $assembly = $organism->get_assembly('2.1_051011');

my $ass_source = $assembly->source;
my $ci = $ass_source->contigs;

my $u_contigs;
my $p_contigs;
my @a_contigs;
my $counter = 0;

my $true = 1;

while($true){
    print $counter."\n";
    $organism->schema->txn_do(
        sub {
            for (1..1000){
                my $contig = $ci->next;
                unless ($contig){
                    $true = 0; 
                    last;
                }
                $counter ++;
                my $id = $contig->id;
                my $name = "Contig".$contig->scaffold->scaffold_num.".".$contig->contig_num;
                my $stated_length = $contig->length;
                my $base_string = $contig->consensus->bases;
                my $string_length = length($base_string);
                my $qualities = $contig->consensus->qualities;
                
                #UPDATING IN DB
                #if ($qualities =~ s/ \*//g){
                #    $contig->consensus->qualities($qualities);
                #    $contig->consensus->update;
                #}
                #if ($stated_length != length($base_string)){
                #    $contig->length($string_length);
                #    $contig->update;
                #}

                #ANALYSIS
                unless ($stated_length and $string_length and $qualities){
                    print "problem with $name, id $id";
                }
                if ($stated_length == $string_length){
                }else{
                print "$stated_length != $string_length\n";
                    $u_contigs++;
                }
                if ($qualities =~ /\*/){
                    print "pads stored!\n" and $p_contigs++;
                }

            }

            print "$u_contigs mismatched length contigs\n" if $u_contigs;
            print "$p_contigs qual strings with pads\n" if $p_contigs;

        }
    );
}

=pod

=head1 NAME
ScriptTemplate - template for new perl script

=head1 SYNOPSIS

=head1 DESCRIPTION 

=cut

#$HeadURL$
#$Id$


