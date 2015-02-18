#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};


use Test::More tests => 3;
use Test::Exception;

use Sub::Install;


use above "Genome";
use Genome::Test::Factory::Build;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Model::ReferenceSequence;
use Genome::Test::Factory::SoftwareResult::User;

my $pkg = 'Genome::Model::Tools::Htseq::Count::Result';

use_ok($pkg);

my $temp_dir = Genome::Sys->create_temp_directory;
my $temp_bam = File::Spec->join($temp_dir, 'raw_all_sequences.bam.sort.bam');
my $temp_sam = File::Spec->join($temp_dir, 'test.sam');
Genome::Sys->write_file(
    $temp_sam,
<<'EOSAM'
@HD	VN:1.0	SO:queryname
@SQ	SN:1	LN:249250621	UR:ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz	AS:GRCh37-lite	M5:1b22b98cdeb4a9304cb5d48026a85128	SP:Homo sapiens
A_READ:1:1:7000:1000	0	1	10000	0	10M999N40M	5	10000	0	AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTT	HHHHHIIIIIGGGGGIIIIIGGGGGIIIIIGGGGGIIIIIGGGGGIIIII
EOSAM
);
Genome::Sys->shellcmd(
    cmd => "samtools view -b -S $temp_sam > $temp_bam",
);

my $refseq = Genome::Test::Factory::Model::ReferenceSequence->setup_reference_sequence_build();
my $anno_model = Genome::Test::Factory::Model::ImportedAnnotation->setup_object();
my $anno_build = Genome::Test::Factory::Build->setup_object(model_id => $anno_model->id, reference_sequence => $refseq);

my @ar;
for(1..3) {
    my $data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(lane => $_);
    $data->sample->extraction_type('rna');
    $data->library->transcript_strand('firststrand');
    my $ar = Genome::InstrumentData::AlignmentResult::PerLaneTophat->__define__(aligner_version => 'test', aligner_params => '', instrument_data => $data, annotation_build => $anno_build);
    $ar->temp_scratch_directory($temp_dir);
    push @ar, $ar;
}

my @params = (
    alignment_results => \@ar,
    app_version => '0.5.4p1',
    result_version => 2,
    limit => 2000,
    minaqual => 1,
    mode => 'intersection-strict',
    users => Genome::Test::Factory::SoftwareResult::User->setup_user_hash(),
);

Sub::Install::reinstall_sub({
    into => 'Genome::Sys',
    as => 'shellcmd',
    code => sub {
        my $class = shift;
        my %params = @_;

        my $cmd = $params{cmd};
        return unless($cmd and $cmd =~ /htseq-count/);

        my $err_file = $params{redirect_stderr};
        return unless $err_file;

        Genome::Sys->write_file(
            $err_file,
            "Things were going so well...\n",
            "An expected error occurred!\n",
        );
        die('Not running htseq-count in this test.');
    },
});

my $htseq_count_method = \&Genome::Model::Tools::Htseq::Count::Result::_run_htseq_count;
my $result;
Sub::Install::reinstall_sub({
    into => $pkg,
    as => '_run_htseq_count',
    code => sub {
        $result = $_[0];
        $result->queue_debug_messages(1);
        $htseq_count_method->(@_);
    },
});

throws_ok(sub { $pkg->get_or_create(@params) }, qr/Not running htseq-count in this test./);
my @debug_messages = $result->debug_messages;

my $last_message = pop @debug_messages;
like($last_message, qr/An expected error occurred!/, 'error log was captured');

