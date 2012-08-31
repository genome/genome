use strict;
use warnings;

use above 'Workflow';

use File::Temp;
use File::Basename;
use Test::More tests => 42;

BEGIN {
    use_ok('GAP::Command');
    use_ok('GAP::Command::GenePredictor::tRNAscan');
}

my $command = GAP::Command::GenePredictor::tRNAscan->create(
                                                            'fasta_file' => File::Basename::dirname(__FILE__).'/data/HPAG1.fasta',
                                                            'domain'     => 'bacteria',
                                                           );

isa_ok($command, 'GAP::Command::GenePredictor');
isa_ok($command, 'GAP::Command::GenePredictor::tRNAscan');

ok($command->execute());

my @features = @{$command->bio_seq_feature()};

ok(@features > 0);

foreach my $feature (@features) {
    isa_ok($feature, 'Bio::SeqFeature::Generic');
}
