#!/usr/bin/env genome-perl

use strict;
use warnings;
#use Genome::Model::phaseSnvs;
use above "Genome"; 
use Test::More tests => 1;
use File::Basename;
#use File::Spec my $dirName = dirname(__FILE__);
use File::Temp qw(tempfile);

#print($INC{"phaseSnvs.pm"}, "\n");

my $plot = Genome::Model::phaseSnvs->create(
                                                         distance  => 700,
                                                         command => "C",
							 vcfFile => "/gscuser/cfederer/Documents/mysnvs.vcf",
                                                         bamFile => "/gscuser/cfederer/Documents/p53.bam",
                                                         chromosome => 17,
							 sample => "none",
							 relax => "false", 
							 snvs => "7579646T7579653C", 
                                                        );
ok($plot->execute,'Executed');
                                                  
