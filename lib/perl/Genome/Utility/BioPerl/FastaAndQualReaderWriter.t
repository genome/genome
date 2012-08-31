#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
require File::Temp;
use Bio::Seq::Quality;

#< Bioseqs >#
my $bioseq = Bio::Seq::Quality->new(
    '-id' => 'bioseq',
    '-seq' => 'AGCT',
    '-qual' => [qw/ 12 13 14 15 /],
);
ok($bioseq, 'Created bioseq');
my $invalid_bioseq_1 = Bio::Seq::Quality->new(
    '-id' => 'invalid_bioseq_1',
    '-seq' => 'AGC1',
    '-qual' => [qw/ 12 13 14 15 /],
);
ok($invalid_bioseq_1, 'Created invalid bioseq');

#< Files >#
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir.');
my $fasta_file = $tmpdir.'/fasta';
my $qual_file = $tmpdir.'/qual';

#< Writer >#
my $writer_class = 'Genome::Utility::BioPerl::FastaAndQualWriter';
use_ok($writer_class) or die;

# Fail - No fasta file
eval{ $writer_class->create; };
ok($@, 'Create  w/o fasta file.');
diag($@);

# Good
my $writer = $writer_class->create(
    fasta_file => $fasta_file,
    qual_file => $qual_file,
);
ok($writer, 'Create for writing.');
ok($writer->write_seq($bioseq), 'Write seq.');
eval{ $writer->write_seq($invalid_bioseq_1); };
ok($@, "Tried to write invalid bioseq: $@");

#< Reader >#
my $reader_class = 'Genome::Utility::BioPerl::FastaAndQualReader';
use_ok($reader_class) or die;

# Fail - No fasta file
eval{ $reader_class->create };
ok($@, 'Create  w/o fasta file.');
diag($@);

my $reader = $reader_class->create(
    fasta_file => $fasta_file,
    qual_file => $qual_file,
);
ok($reader, 'Create for reading.');
my $same_bioseq = $reader->next_seq;
ok($same_bioseq, 'Got bioseq from reader.');
is($bioseq->id, $same_bioseq->id, 'Compared bioseqs ids from reader and writer');
is($bioseq->seq, $same_bioseq->seq, 'Compared bioseqs seqs from reader and writer');

done_testing();


=pod

=head1 Tests

=head1 Disclaimer

 Copyright (C) 2010 Washington University Genome Sequencing Center

 This script is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
