#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::Trim::ByQual') or die;

# Create failures
ok(scalar(Genome::Model::Tools::Sx::Trim::ByQual->create->__errors__), 'Create fails w/o quality');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByQual->create(quality => 'all')->__errors__), 'Create fails w/ quality => all');
ok(scalar(Genome::Model::Tools::Sx::Trim::ByQual->create(quality => -1)->__errors__), 'Create fails w/ quality => -1');

# Trimming
my $trimmer = Genome::Model::Tools::Sx::Trim::ByQual->create(quality => 2);
ok($trimmer, 'Create trimmer');
my ($seq, $qual);
for (0..93) { $seq .= 'A'; $qual .= chr($_ + 33); }
$qual = reverse $qual;
is(length($seq), length($qual), 'Seq and qual length match');
my $seq1a = { seq => $seq, qual => $qual };
my $seq1b = { seq => $seq, qual => $qual };
my $evaluator = $trimmer->_create_evaluator;
$evaluator->([ $seq1a, $seq1b ]);
is(length($seq1a->{seq}), 91, 'Seq 1a seq trimmed to 91 bases');
is(length($seq1a->{qual}), 91, 'Seq 1a qual trimmed to 91 bases');
is(length($seq1b->{seq}), 91, 'Seq 1b seq trimmed to 91 bases');
is(length($seq1b->{qual}), 91, 'Seq 1b qual trimmed to 91 bases');
my $seq2 = { seq => $seq, qual => $qual };
$trimmer->quality(25);
$evaluator = $trimmer->_create_evaluator;
$evaluator->([ $seq2 ]);
is(length($seq2->{seq}), 68, 'Seq 2 seq trimmed to 68 bases');
is(length($seq2->{qual}), 68, 'Seq 2 qual trimmed to 68 bases');
$trimmer->quality(100);
$evaluator = $trimmer->_create_evaluator;
$evaluator->([ $seq2 ]);
is(length($seq2->{seq}), 0, 'Seq 3 seq trimmed to 0 bases');
is(length($seq2->{qual}), 0, 'Seq 3 qual trimmed to 0 bases');

done_testing();
exit;

