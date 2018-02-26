#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More skip_all => 'broken, but this is not in production yet, so it is disabled';

use_ok('Genome::WorkOrderItem');

#
# Caveats
#  USING REAL PROD DATA!
#  NOT testing a work order that does not have seq products
#   Couldn't find one easily
#  Numbers of inst data models can change, so not testing exact count or names/ids
#   It is also possible that hesee model may get deleted, in which case this test will
#   fail.
#

#< This is a work order item that has 454 data runs and a model
# woi_id
#  141
# inst data ids 
#  2737764362 RunRegion454
#  2737764363 RunRegion454
#  2852539831 RegionIndex454
#  2852539832 RegionIndex454
# model
#  2744704120
my $woi = Genome::WorkOrderItem->get(141);
ok($woi, 'Got work order item (141) w/ 454 indexed regions');

my @sequence_products = $woi->sequence_products;
ok(@sequence_products, 'sequence products');
#print Dumper(@sequence_products);

# FIXME The instrument_data method works, but not for index region 454.  
#       But I can get these inst data when run on the command line:
#       perl -M'above "Genome"' -e 'while (<>) { chomp; my $i = Genome::InstrumentData->get(2737764362); print("$_ => ", ref($i) || "none", "\n"); }'
#
#       Not sure what the problem is
#
#my @instrument_data = $woi->instrument_data;
#cmp_ok(@instrument_data, '>=', 4, 'instrument data'); # >= incase more is added
#print Dumper(@instrument_data);

my @models = $woi->models;
ok(@models, 'models for 454 work order');
#print Dumper(\@models);
#>

#< This is a work order item that has 8 sanger runs and a model.  Testing
#   sanger in addition to 454 because woi track reads, not runs, and we
#   assign runs to models. So this gets the reads, then the runs.  Where
#   454/solexa is stored directly.
# woi id
#  2642
# inst data ids
#  21jul09.909pmcb1
#  21jul09.906pmab2
#  21jul09.906pmaa2
#  21jul09.906pmaa1
#  21jul09.909pmca1
#  21jul09.909pmca2
#  21jul09.909pmcb2
#  21jul09.906pmab1
# model
#   2852891065
$woi = Genome::WorkOrderItem->get(2642);
ok($woi, 'Got work order item (2642) w/ sanger reads');

@sequence_products = $woi->sequence_products;
cmp_ok(@sequence_products, '>=', 8, 'got woi sequence products'); # >= incase more is added
#print Dumper({ map { $_->prep_group_id => 1 } @sequence_products }); # reads

# NOT TESTING instrument_data here.  The 454 test above is sufficient.

@models = $woi->models;
ok(@models, 'models for sanger work order');
#print Dumper(\@models);
#>

done_testing();


=pod

=head1 Tests

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This script is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut
