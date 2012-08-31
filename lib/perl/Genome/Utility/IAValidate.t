
use Test::More tests => 8;
use above 'Genome';

use File::Slurp;

BEGIN {
    use_ok("Genome::Utility::IAValidate");
}


my $prot_line1 = "56107\t66199\tENSP00000374929\tMWGAFLLYVSMKMGGTAGQSLEQPSEVTAVEGAIVQINCTYQTSGFYGLSWYQQHDGGAPTFLSYNALDGLEETGRFSSFLSRSDSYGYLLLQELQMKDSASYFCAVR";

my $prot_line2 = "12345\t54321\tENSPTHISISTEST\t";

ok(Genome::Utility::IAValidate::check_line($prot_line1,'protein'));
my $retstring = Genome::Utility::IAValidate::check_line($prot_line2,'protein');
ok($retstring eq
   "12345\t54321\tENSPTHISISTEST\t:no peptide sequence, ENSPTHISISTEST");


# hit transcripts, genes?
my $transcript1 = "117851\t67956\t1873\t3533\tENST00000404059\tensembl\tknown\t+1\t1";
my $transcript2 = "117851\t67956\t1\t1\tENST00000404059\tensembl\tknown\t+1\t1";
ok(Genome::Utility::IAValidate::check_line($transcript1,'transcript'));
my $retstring_transcript = Genome::Utility::IAValidate::check_line($transcript2,
   'transcript');
#print $retstring_transcript,"\n";
ok($retstring_transcript eq "117851\t67956\t1\t1\tENST00000404059\tensembl\tknown\t+1\t1:start/stop same,");


# hit transcript sub structures

my $tss1 = "736555\t66206\tcds_exon\t21181033\t21181075\t1\t1\tATGTGGGGAGTTTTCCTTCTTTATGTTTCCATGAAGATGGGAG";
my $tss2 = "736553\t66206\tflank\t21131033\t21181032\t1\t\t";
my $tss3 = "739290\t58931\tflank\t-48396\t1603\t1\t\t";
ok(Genome::Utility::IAValidate::check_line($tss1,'sub_structure'));
ok(Genome::Utility::IAValidate::check_line($tss2,'sub_structure'));
ok(Genome::Utility::IAValidate::check_line($tss3,'sub_structure') eq
   "739290\t58931\tflank\t-48396\t1603\t1\t\t:start position is negative");
# _fin_

# $Id$
