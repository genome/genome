use strict;
use warnings;

use above "PAP";
use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use File::Temp;
use File::Basename;
use Test::More tests => 532;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::BlastP');
}

my $tempdir = File::Temp::tempdir(
                                  'PAP_BlastP_test_XXXXXXXX',
                                  DIR     => '/tmp',
                                  CLEANUP => 1,
                                 );

my $command = PAP::Command::BlastP->create(
                                           'fasta_file'      => File::Basename::dirname(__FILE__).'/data/B_coprocola.chunk.fasta', 
                                           'report_save_dir' => $tempdir, 
                                          );
                                          
isa_ok($command, 'PAP::Command::BlastP');

SKIP: {
    skip "long test; run manually setting RUNBLASTP=1", 529 unless $ENV{RUNBLASTP};
ok($command->execute());

my $ref = $command->bio_seq_feature();

is(ref($ref), 'ARRAY');

foreach my $feature (@{$ref}) {

    isa_ok($feature, 'Bio::SeqFeature::Generic');

    my $annotation_collection = $feature->annotation();

    isa_ok($annotation_collection, 'Bio::Annotation::Collection');

    my ($annotation) = $annotation_collection->get_Annotations();

    if (defined($annotation)) {
    
        isa_ok($annotation, 'Bio::Annotation::DBLink');
        is($annotation->database(), 'GenBank');
        like($annotation->primary_id(), qr/^\w+$/);

        ok($feature->has_tag('blastp_bit_score'));
        ok($feature->has_tag('blastp_evalue'));
        ok($feature->has_tag('blastp_percent_identical'));
        ok($feature->has_tag('blastp_query_start'));
        ok($feature->has_tag('blastp_query_end'));
        ok($feature->has_tag('blastp_subject_start'));
        ok($feature->has_tag('blastp_subject_end'));
        ok($feature->has_tag('blastp_hit_name'));
        ok($feature->has_tag('blastp_hit_description'));

        my ($hit_description) = $feature->each_tag_value('blastp_hit_description');
        
        unlike($hit_description, qr/>/);
    
    }
    
    ok($feature->has_tag('blastp_category'));

    if ($feature->display_name() eq 'BACCOPFNL_Contig86.GeneMark.10') {
        my ($blastp_category) = $feature->each_tag_value('blastp_category');    
        is($blastp_category, 'Hypothetical Protein', 'blastp_category is plain hypothetical');
    }
    elsif ($feature->display_name() eq 'BACCOPFNL_Contig86.GeneMark.14') {
         my ($blastp_category) = $feature->each_tag_value('blastp_category');
         like($blastp_category, qr/^Hypothetical Protein similar to/, 'blastp_category is hypothetical with similarity');
    }

    my $display_name = $feature->display_name();

    ok(-e "$tempdir/$display_name.blastp.bz2", "archived blast report exists for $display_name");
    is(system("bzcat /$tempdir/$display_name.blastp.bz2"), 0, "bzcat can read archived blast report for $display_name");
    
}

my $blast_report_fh = $command->blast_report();

isa_ok($blast_report_fh, 'File::Temp');
ok($blast_report_fh->opened());

}

