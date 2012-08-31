use strict;
use warnings;

use above "PAP";
use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use File::Temp;
use Test::More tests => 4;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::BerBlastp');
}

my @seqnames = ('SHORTTESTDFT2_Contig1.1.GeneMark.p5_hybrid.10', 
                'SHORTTESTDFT2_Contig1.1.GeneMark.p5_hybrid.11',
                'SHORTTESTDFT2_Contig1.1.GeneMark.p5_hybrid.12',) ;
my $blastpquery = "/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/data/panda/AllGroup/AllGroup.niaa";
my $fastadir = "/gsc/var/cache/testsuite/data/PAP-Command-BerBlastp/fasta";
my $berdirpath = File::Temp::tempdir("ber-blastp-berdir-XXXXXX", CLEANUP => 1);
#my $berdirpath = File::Temp::tempdir("/tmp/ber-blastp-berdir-XXXXXX");

my $b = PAP::Command::BerBlastp->create(
                      sequence_names => \@seqnames,
                      fastadir       => $fastadir,
                      blastp_query   => $blastpquery,
                      berdirpath     => $berdirpath,
                     );

isa_ok($b,'PAP::Command::BerBlastp');

ok($b->execute(),'we can run the blastp job(s)');
# there is'nt too much to test here - parsing or conversion doesn't take
# place in this module.
