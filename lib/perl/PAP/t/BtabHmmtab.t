use strict;
use warnings;
use above 'PAP';

use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use File::Temp qw(tempdir);
use File::Find;
use Test::More tests => 26;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::BtabHmmtab');
}

my @seqs = ( );
our @files;
my $locus_tag = 'SHORTTESTDFT2';

# original location
my $ofastadir = '/gsc/var/cache/testsuite/data/PAP-Command-BtabHmmtab/btabhtab-test/fasta';
my $oberdirpath = '/gsc/var/cache/testsuite/data/PAP-Command-BtabHmmtab/btabhtab-test/ber';
my $ohmmdirpath = '/gsc/var/cache/testsuite/data/PAP-Command-BtabHmmtab/btabhtab-test/hmm';
my $srcdirpath = '/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/src/';

my $dir = tempdir('/tmp/pap-test-XXXXXXX', CLEANUP => 1);
#my $dir = tempdir('/tmp/pap-test-XXXXXXX');
diag($dir);
my $fastadir = $dir."/btabhtab-test/fasta";
my $berdirpath = $dir."/btabhtab-test/ber";
my $hmmdirpath = $dir."/btabhtab-test/hmm";
my $bsubfilepath = $dir.'/btabhtab-test/bsubERRfiles';
unless(-d $dir) {
    mkdir($dir);
    mkdir($fastadir);
    mkdir($berdirpath);
    mkdir($hmmdirpath);
    mkdir($bsubfilepath);
}
my $rv = system("cp -pr /gsc/var/cache/testsuite/data/PAP-Command-BtabHmmtab/btabhtab-test/ $dir");




my $b = PAP::Command::BtabHmmtab->create(
        locus_tag => $locus_tag,
        fastadir => $fastadir,
        berdirpath => $berdirpath,
        hmmdirpath => $hmmdirpath,
        srcdirpath => $srcdirpath,
        bsubfiledirpath => $bsubfilepath,
        sequence_names => \@seqs, # is this needed????
);

isa_ok($b,'PAP::Command::BtabHmmtab');

my @fasta_store = $b->get_fasta_store;
is(scalar(@fasta_store), 9, 'right number of fasta files(get_fasta_store)');
my @ber_store = $b->get_ber_store;
is(scalar(@ber_store), 9, 'right number of blastp results(get_ber_store)');
my @hmm_store = $b->get_hmm_store;
is(scalar(@hmm_store), 9, 'right number of hmmpfam results(get_hmm_store)');

ok($b->ber2btab(),'ber2btab runs');
ok($b->hmm2htab(),'hmm2htab runs');

# need to check the output files...
# eventually have to get rid of output files.
find(\&wanted, $hmmdirpath, $berdirpath);

foreach my $file (@files)
{
    ok(-f $file, 'result file exists');
    # md5sum of output?
    unlink($file);
}



sub wanted 
{
    if ($_ =~ /(\.btab|\.htab)$/) {
        push(@files, $File::Find::name);
    }
}
