#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More;

#plan skip_all => 'This test uses the real instrument solexa data, which is not desirable. But this need to be run to test new class if any';

my $id = 2889838607;
my $old = Genome::Site::TGI::InstrumentData::Solexa->get($id);
ok($old, 'got old index illumina');
my $new = Genome::Site::TGI::Synchronize::Classes::IndexIllumina->get($id);
ok($new, 'got index illumina');
my $old_data = $old->{db_committed};
my $new_data = $new->{db_committed};
for my $prop (qw/ sequencing_platform subclass_name project_name /) {
    delete $old_data->{$prop};
}
is_deeply($new_data, $old_data, 'new and old match');
#print Data::Dumper::Dumper($new_data, $old_data);

my ($old_dp, $old_ip) = Genome::Site::TGI::Synchronize::UpdateApipeClasses->_get_direct_and_indirect_properties_for_object(
    $old,
    'Genome::InstrumentData::Solexa',
    qw/ sample_name sample_id /
);
my ($new_dp, $new_ip) = Genome::Site::TGI::Synchronize::UpdateApipeClasses->_get_direct_and_indirect_properties_for_object(
    $new,
    'Genome::InstrumentData::Solexa',
    qw/ sample_name sample_id /
);
is_deeply($new_dp, $old_dp, 'new and old direct properties match');
is_deeply($new_ip, $old_ip, 'new and old direct properties match');
#print Dumper($new_dp, $old_dp, $new_ip, $old_ip);

done_testing();
