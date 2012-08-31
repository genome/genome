use strict;
use warnings;

use above 'PAP';
use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use Cwd;
use File::Temp qw/ tempdir /;
use Test::More tests => 8;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::Hmmpfam');
}

my @sequence_names = (
        'SHORTTESTDFT2_Contig108.1.GeneMark.p5_hybrid.62',
        'SHORTTESTDFT2_Contig108.1.GeneMark.p5_hybrid.64',
        'SHORTTESTDFT2_Contig108.1.GeneMark.p5_hybrid.9',
        'SHORTTESTDFT2_Contig108.1.Glimmer3.p5_hybrid.18',
);
my $fastadir = '/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/data/genomes/SHORTTESTDFT2/fasta/';
my $hmmdatabase = "/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/data/ALL_LIB.HMM";

SKIP: {
    skip "test data missing? test data needs to be put somewhere better", 6 unless 0;
my $hmmdir = tempdir( CLEANUP => 1);
my $workdir = $hmmdir;
diag($hmmdir);
my $h = PAP::Command::Hmmpfam->create(
        hmm_database   => $hmmdatabase,
        hmmdirpath     => $hmmdir,
        fasta_dir      => $fastadir,
        sequence_names => \@sequence_names,
        workdir        => $workdir,
);

isa_ok($h, 'PAP::Command::Hmmpfam');

diag('execution takes a few minutes');
ok($h->execute(), 'running hmmpfam');

foreach my $seqname (@sequence_names)
{
    my $hmmpfamout = $hmmdir."/".$seqname.".hmmpfam";
    ok(-f $hmmpfamout, 'output file '.$hmmpfamout.' exists');
}
}
