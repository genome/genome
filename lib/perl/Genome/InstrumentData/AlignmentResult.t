#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 6;

# This test was auto-generated because './InstrumentData/AlignmentResult.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
my $class = 'Genome::InstrumentData::AlignmentResult';
use_ok($class);

class Genome::InstrumentData::AlignmentResultTester {
    is => $class,
};
my $inst_data = Genome::InstrumentData::Imported->__define__(id => -111);
my $alignment_result = Genome::InstrumentData::AlignmentResultTester->__define__(
    id => -1337,
    output_dir => '/dev/null',
    instrument_data_id => $inst_data->id,
    aligner_name => 'bwa',
    aligner_version => '1',
    aligner_params => 'NA',
);
ok($alignment_result, 'defined alignment result');

#  define qc result
my $qc_result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->__define__(alignment_result_id => $alignment_result->id);
ok($qc_result, 'define qc result for alienment result');
ok($alignment_result->add_user(user => $qc_result, label => 'uses'), 'add qc result as user of alignment result');

# delete
ok($alignment_result->delete, );
ok(ref($qc_result) eq 'UR::DeletedRef', 'deleted qc result');

done_testing();
