#!/usr/bin/env genome-perl

use strict;
use warnings;

use Data::Dumper;
use IO::File; 

my $f1 = IO::File->new("< Schema.new.pm") or die "$!\n";
unlink "Schema.eb.pm";
my $f2 = IO::File->new("> Schema.eb.pm") or die "$!\n"; 
while ( my $line = $f1->getline )
{
    # $assembly_cache->set('assembly', 'line', $as_line);
    if ($line =~ /cache\-\>set\('/)
    {
        $line =~ s/set\(\'/\{/; 
        $line =~ s/', '/\}\-\>\{/;
        $line =~ s/', /} = /; 
        $line =~ s/\)//; 
    } 
    # $contig_cache->set($contig_name, 'bs_start', $start);
    elsif ($line =~ /cache\-\>set\(\$/)
    {
        $line =~ s/set\(/\{/; 
        $line =~ s/, '/\}\-\>\{/;
        $line =~ s/', '?/} = /; 
        $line =~ s/'?\)//; 
    } 
    # $assembly_cache->get('assembly', 'line');
    # $contig_cache->get($contig_name, 'bs_start');
    elsif ($line =~ /cache\-\>get\('/)
    {
        $line =~ s/get\(\'/\{/; 
        $line =~ s/', '/\}\-\>\{/;
        $line =~ s/', /} = /; 
        $line =~ s/\)//; 
    } 
    # $contig_cache->set($contig_name, 'bs_start', $start);
    elsif ($line =~ /cache\-\>get\(\$/)
    {
        $line =~ s/get\(/\{/; 
        $line =~ s/, '/\}\-\>\{/;
        $line =~ s/', '?/} = /; 
        $line =~ s/'?\)//; 
    } 

    $f2->print($line);
}

$f1->close;
$f2->close;

=pod

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

