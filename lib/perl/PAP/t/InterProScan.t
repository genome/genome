use strict;
use warnings;

use above "PAP";
use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use File::Temp;
use File::Basename;
#use Test::More tests => 289;
use Test::More;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::InterProScan');
}

my $tempdir = File::Temp::tempdir(
		'PAP_InterPro_test_XXXXXXXX',
		DIR     => '/tmp',
		CLEANUP => 1,
		);

my $command = PAP::Command::InterProScan->create(
		'fasta_file'      => File::Basename::dirname(__FILE__).'/data/B_coprocola.chunk.fasta',
		'report_save_dir' => $tempdir,
		'locus_tag'	=> 'TEST',
		);
isa_ok($command, 'PAP::Command::InterProScan');

ok($command->execute());

my $ref = $command->bio_seq_feature();

is(ref($ref), 'ARRAY');

foreach my $feature (@{$ref}) {

    isa_ok($feature, 'Bio::SeqFeature::Generic');

    ok($feature->has_tag('interpro_analysis'));
    ok($feature->has_tag('interpro_evalue'));
    ok($feature->has_tag('interpro_description'));

    my $annotation_collection = $feature->annotation();

    isa_ok($annotation_collection, 'Bio::Annotation::Collection');

    my @annotations = $annotation_collection->get_Annotations();

    foreach my $annotation (@annotations) {

        isa_ok($annotation, 'Bio::Annotation::DBLink'); 

    }

}

ok(-e "$tempdir/merged.raw.sorted.bz2", "archive interpro sorted output exists");
ok(! -z "$tempdir/merged.raw.sorted.bz2", "archive interpro sorted output is not empty");
ok(-e "$tempdir/merged.raw.bz2", "archive interpro raw output exists");
ok(! -z "$tempdir/merged.raw.bz2", "archive interpro raw output is not empty");
is(system("bzcat $tempdir/merged.raw.sorted.bz2 > /dev/null"), 0, "bzcat can read archived sorted output");
is(system("bzcat $tempdir/merged.raw.bz2 > /dev/null"), 0, "bzcat can read archived raw output");

my $interpro_output_fh = $command->iprscan_output();

isa_ok($interpro_output_fh, 'File::Temp');
ok($interpro_output_fh->opened());
done_testing();
