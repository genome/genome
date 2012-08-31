#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 63;

BEGIN {	
    use_ok('Genome::Model::Tools::RefCov::ROI::RegionI');
    use_ok('Genome::Model::Tools::RefCov::ROI::Region');
}

my $region1 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 1,
    end => 1000,
    strand => 1
);

my $region2 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 100,
    end => 900,
    strand => -1
);

my $subtracted = $region1->subtract($region2);
ok(defined($subtracted));
is(scalar(@$subtracted), 2);
foreach my $region (@$subtracted) {
    ok($region->start == 1 || $region->start == 901);
    ok($region->end == 99 || $region->end == 1000);
}

$subtracted = $region2->subtract($region1);
ok(!defined($subtracted));
$subtracted = $region1->subtract($region2, 'weak');
ok(!defined($subtracted));
$subtracted = $region1->subtract($region2, 'strong');
ok(!defined($subtracted));

my $region3 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 500,
    end => 1500,
    strand => 1
);

$subtracted = $region1->subtract($region3);
ok(defined($subtracted));
is scalar(@$subtracted), 1;
my $subtracted_i = @$subtracted[0];
is($subtracted_i->start, 1);
is($subtracted_i->end, 499);


my $genome_region = Genome::Model::Tools::RefCov::ROI::Region->create(
    start =>10,
    end =>20,
    strand =>1
);

isa_ok($genome_region,'Genome::Model::Tools::RefCov::ROI::Region', 'Genome::Model::Tools::RefCovRegion object');
is($genome_region->strand, 1);

my $genome_region2 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start=>15,
    end=>25,
    strand=>1
);

isa_ok($genome_region2,'Genome::Model::Tools::RefCov::ROI::Region', 'Genome::Model::Tools::RefCovRegion object');
is($genome_region2->strand, 1);

my $r = Genome::Model::Tools::RefCov::ROI::Region->create();
is ( $r->strand(0), 0 ) ;
is ( $r->start(27), 27 );
is ( $r->end(28), 28 ) ;

ok(! defined $r->intersection($genome_region2));

$r = $genome_region->union($genome_region2);
is($r->start, 10);
is($r->end, 25);

$r = $genome_region->intersection($genome_region2);
is ( $r->start, 15  ) ;
is ( $r->end, 20    );
is ( $r->strand, 1  );

# intersection and union can also take lists
my $genome_region3 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start=>18,
    end=>30
);
isa_ok($genome_region3,'Genome::Model::Tools::RefCov::ROI::Region', 'Genome::Model::Tools::RefCovRegion object');

$r = $genome_region->intersection([$genome_region2, $genome_region3]);
ok( ( $r->start == 18 ) && ( $r->end == 20 ));
$r = Genome::Model::Tools::RefCov::ROI::Region->intersection([$genome_region, $genome_region2, $genome_region3]);
ok($r->start == 18 && $r->end == 20);
$r = $genome_region->union($genome_region2, $genome_region3);
ok( ( $r->start == 10 ) && ( $r->end == 30 ) );
$r = Genome::Model::Tools::RefCov::ROI::Region->union($genome_region, $genome_region2, $genome_region3);
ok( ( $r->start == 10 ) && ( $r->end == 30 ) );
$genome_region3->start(21);
ok (! $genome_region->intersection([$genome_region2, $genome_region3]));

ok (! $genome_region->contains($genome_region2));
ok (! $genome_region2->contains($genome_region));
ok ($genome_region->overlaps($genome_region2));
ok ($genome_region2->overlaps($genome_region));

# testing strand
$genome_region3 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 15,
    end => 25,
    strand => 1
);

my $genome_region4 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 15,
    end => 25,
    strand => -1
);

isa_ok($genome_region4,'Genome::Model::Tools::RefCov::ROI::Region', 'Genome::Model::Tools::RefCovRegion object');

my $genome_region5 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 15,
    end => 25,
    strand => 0
);

isa_ok($genome_region5,'Genome::Model::Tools::RefCov::ROI::Region', 'Genome::Model::Tools::RefCovRegion object');

my $genome_region6 = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 20,
    end => 30,
    strand => -1
);

isa_ok($genome_region6,'Genome::Model::Tools::RefCov::ROI::Region', 'Genome::Model::Tools::RefCovRegion object');

ok $genome_region3->_ignore($genome_region4), ' 1 & -1' ;     
ok $genome_region3->_weak($genome_region3),' 1 & 1 true' ;       
ok $genome_region3->_weak($genome_region5), ' 1 & 0 true' ;       
ok (! $genome_region3->_weak($genome_region4), ' 1 & -1 false' );   
ok $genome_region3->_strong($genome_region3), ' 1 & 1 true' ;     
ok (! $genome_region3->_strong($genome_region5), ' 1 & 0 false' ); 
ok (! $genome_region3->_strong($genome_region4), ' 1 & -1 false' ); 

ok ! ( $genome_region3->overlaps($genome_region4,'weak'));
ok ! ( $genome_region4->overlaps($genome_region3,'weak'));
ok ! ( $genome_region3->overlaps($genome_region4,'strong')); 
ok ! ( $genome_region4->overlaps($genome_region3,'strong')); 

$genome_region3->strand(0);

ok  ( $genome_region3->overlaps($genome_region4,'weak'));
ok  ( $genome_region4->overlaps($genome_region3,'weak')); 
ok ! ( $genome_region3->overlaps($genome_region4,'strong'));
ok ! ( $genome_region4->overlaps($genome_region3,'strong')); 

# if strands are different then intersection() should return 0...
$r = $genome_region3->intersection($genome_region4);
is ( $r->strand, 0 );

# or if both strands are -1 then -1 should be returned
$r = $genome_region6->intersection($genome_region4);
is ( $r->strand, -1 );

# test implemention of offsetStranded:
$r = Genome::Model::Tools::RefCov::ROI::Region->create(
    start => 30,
    end => 40,
    strand => -1
);
isa_ok($r, 'Genome::Model::Tools::RefCov::ROI::Region', 'Genome::Model::Tools::RefCov::ROI::Region object') ;
is ($r->offsetStranded(-5,10)->toString, '(20, 45) strand=-1');
is ($r->offsetStranded(+5,-10)->toString, '(30, 40) strand=-1');
$r->strand(1);
is ($r->offsetStranded(-5,10)->toString, '(25, 50) strand=1');
is ($r->offsetStranded(+5,-10)->toString, '(30, 40) strand=1');
