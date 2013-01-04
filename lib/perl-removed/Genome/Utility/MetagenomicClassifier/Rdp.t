#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Bio::SeqIO;
use Test::More;

#this serves as the test case for all versions of Rdp - inheriting from Genome::Utility::MetagenomicClassifier::Rdp

use_ok('Genome::Utility::MetagenomicClassifier::Rdp::Version2x1');
use_ok('Genome::Utility::MetagenomicClassifier::Rdp::Version2x2');

#each Version has to be run separately to avoid running the same JVM instance

done_testing();


#$HeadURL$
#$Id$
