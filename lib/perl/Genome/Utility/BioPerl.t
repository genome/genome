#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Bio::Seq;
require Bio::Seq::Quality;
require File::Temp;
use Test::More;

my $class = 'Genome::Utility::BioPerl';
use_ok($class) or die;

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'Created temp dir.');
my $fasta_file = $tmpdir.'/fasta';

#< Bioseq Validation >#
# ok - Bio::Seq
my $bioseq = Bio::Seq->new(
    '-id' => 'bioseq',
    '-seq' => 'AGCT',
);
ok($bioseq, 'Created bioseq');
ok($class->validate_bioseq($bioseq), 'validate Bio::Seq');
# ok - Bio::Seq::Quality
$bioseq = Bio::Seq::Quality->new(
    '-id' => 'bioseq',
    '-seq' => 'AGCT',
    '-qual' => [qw/ 12 13 14 15 /],
);
ok($bioseq, 'Created bioseq');
ok($class->validate_bioseq($bioseq), 'validate Bio::Seq::Quality');
# fail - bioseq undef
eval{ $class->validate_bioseq(); }; 
ok($@, 'validate fail - bioseq undef');
diag($@);
# fail - not a bioseq
eval{ $class->validate_bioseq(1); }; 
ok($@, 'validate fail - not a Bio::Seq');
diag($@);
# fail - bad chars
$bioseq = Bio::Seq->new(
    '-id' => 'bioseq',
    '-seq' => 'AGCT==',
);
ok($bioseq, 'Created bioseq');
eval{ $class->validate_bioseq($bioseq); }; 
ok($@, 'validate fail - bad chars');
diag($@);
# fail - uneven length of seq and qual
$bioseq = Bio::Seq::Quality->new(
    '-id' => 'bioseq',
    '-seq' => 'AGCT',
    '-qual' => [qw/ 12 13 14 15 1 1 1 1 1 /],
);
ok($bioseq, 'Created bioseq');
eval{ $class->validate_bioseq($bioseq); };
ok($@, 'validate fail - uneven length of seq and qual');
diag($@);

# readers/writers
eval{ $class->create_bioseq_writer(undef); };
ok($@, 'can\'t create bioseq writer w/o file');
diag($@);
eval{ $class->create_bioseq_reader(undef); };
ok($@, 'can\'t create bioseq reader w/o file');
diag($@);
eval{ $class->create_bioseq_reader($fasta_file); };
ok($@, 'can\'t create bioseq reader w/o existing file');
diag($@);
ok($class->create_bioseq_writer($fasta_file), 'created bioseq writer');
ok($class->create_bioseq_reader($fasta_file), 'created bioseq reader');

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
