use strict;
use warnings;

use above 'GAP';

use File::Temp;
use File::Basename;
use Test::More tests => 51;

BEGIN {
    use_ok('GAP::Command');
    use_ok('GAP::Command::GenePredictor::RfamScan');
}

my $command = GAP::Command::GenePredictor::RfamScan->create(
                                                            'fasta_file' => File::Basename::dirname(__FILE__).'/data/HPAG1.fasta',
                                                           );

isa_ok($command, 'GAP::Command::GenePredictor');
isa_ok($command, 'GAP::Command::GenePredictor::RfamScan');

SKIP: {
    skip "takes a long long time; run manually with RUNRFAM=1", 47 unless $ENV{RUNRFAM};
ok($command->execute());

my @features = @{$command->bio_seq_feature()};

ok(@features > 0);

foreach my $feature (@features) {
    isa_ok($feature, 'Bio::SeqFeature::Generic');
}
}
