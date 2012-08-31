use strict;
use warnings;

use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use Cwd;
use File::Temp;
use File::Basename;
#use Test::More tests => 81;
use Test::More qw(no_plan);

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::KEGGScan');
}

my $tempdir = File::Temp::tempdir(
                                  'PAP_KEGGscan_test_XXXXXXXX',
                                  DIR     => '/tmp',
                                  CLEANUP => 1,
                                 );

my $command = PAP::Command::KEGGScan->create(
                                             'fasta_file'      => File::Basename::dirname(__FILE__).'/data/B_coprocola.chunk.fasta',
                                             'report_save_dir' => $tempdir,
                                            );
isa_ok($command, 'PAP::Command::KEGGScan');

ok($command->execute());

my $ref = $command->bio_seq_feature();

is(ref($ref), 'ARRAY');

foreach my $feature (@{$ref}) {

    isa_ok($feature, 'Bio::SeqFeature::Generic');
     
    ok($feature->has_tag('kegg_evalue'));
    ok($feature->has_tag('kegg_description'));

    my ($evalue) = $feature->each_tag_value('kegg_evalue');
    
    cmp_ok($evalue, '<', '0.01');
    
    my $annotation_collection = $feature->annotation();

    isa_ok($annotation_collection, 'Bio::Annotation::Collection');

    my @annotations = $annotation_collection->get_Annotations();
    my @dblinks     = grep { $_->isa('Bio::Annotation::DBLink') } @annotations;
    
    my ($gene_dblink)      = grep { $_->primary_id() =~ /^\w{3}\:\w+$/ } @dblinks;
    my ($orthology_dblink) = grep { $_->primary_id() =~ /^K\d+$/       } @dblinks;
  
}

ok(-e "$tempdir/REPORT-top.ks.bz2", "archive KEGGScan top output exists");
ok(! -z "$tempdir/REPORT-top.ks.bz2", "archive KEGGScan top output is not zero byte");
ok(-e "$tempdir/REPORT-full.ks.bz2", "archive KEGGScan full output exists");
ok(! -z "$tempdir/REPORT-full.ks.bz2", "archive KEGGScan full output is not zero byte");
is(system("bzcat $tempdir/REPORT-top.ks.bz2 > /dev/null"), 0, 'bzcat can read archived raw top output');
is(system("bzcat $tempdir/REPORT-full.ks.bz2 > /dev/null"), 0, 'bzcat can read archived raw full output');

#done_testing();
