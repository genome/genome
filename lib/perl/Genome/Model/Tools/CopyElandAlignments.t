#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Test::More tests => 2;

# AUTHOR : Todd Wylie  <twylie@wustl.edu>
# DATE   : Tue Mar 30 10:41:53 CDT 2010

# TEST (1)
use_ok( 'Genome::Model::Tools::CopyElandAlignments' );

# TEST (2)
my $myTest = Genome::Model::Tools::CopyElandAlignments->create(
                                                               flowcell_id => '12345',
                                                               outdir      => '/tmp/impossible_to_find_directory',
                                                               fastq       => 1,
                                                               lanes       => [1, 2, 3],
                                                              );
isa_ok(
       $myTest,
       'Genome::Model::Tools::CopyElandAlignments'
      );

__END__
