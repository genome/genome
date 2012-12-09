#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Data::Dumper 'Dumper';
use File::Grep 'fgrep';
require File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::MetagenomicClassifier::Rdp') or die;

# tmp output dir and file
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
ok(-d $tmpdir, 'created tmp dir');
my $tmp_rdp_file = $tmpdir.'/U_PR-JP_TS1_2PCA.fasta.rdp';
my $fasta = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-MetagenomicClassifier/U_PR-JP_TS1_2PCA.fasta';
ok(-s $fasta, 'Test fasta exists');

# create and execute
my $rdp = Genome::Model::Tools::MetagenomicClassifier::Rdp->create(
        input_file => $fasta,
        output_file => $tmp_rdp_file,
        training_set => 'broad',
        version => '2x1',
        metrics => 1,
);
ok($rdp, 'Created rdp classifier');
$rdp->dump_status_messages(1);
ok($rdp->execute, 'Execute rdp classifier');

# compare output
my $fa_reader = Genome::Model::Tools::Sx::PhredReader->create(
    file => $fasta,
);
ok($fa_reader, 'create fasta reader');
my $cl_reader = Genome::Model::Tools::MetagenomicClassifier::ClassificationReader->create(
    file => $tmp_rdp_file,
);
ok($cl_reader, 'create reader') or die;
my @classifications;
while ( my $fasta = $fa_reader->read ) {
    my $classification = $cl_reader->read;
    is($fasta->{id}, $classification->{id}, 'id matches');
    push @classifications, join('-', map { $classification->{$_}->{id} } (qw/ domain phylum order class /));
}
is_deeply(
    \@classifications, 
    [qw/
        Bacteria-Bacteroidetes-Bacteroidales-Bacteroidetes
        Bacteria-Firmicutes-Clostridiales-Clostridia
        Bacteria-Firmicutes-Clostridiales-Clostridia
        Bacteria-Verrucomicrobia-Verrucomicrobiales-Verrucomicrobiae
        Bacteria-Firmicutes-Clostridiales-Clostridia
        Bacteria-Firmicutes-Clostridiales-Clostridia
        Bacteria-Verrucomicrobia-Verrucomicrobiales-Verrucomicrobiae
        Bacteria-Firmicutes-Clostridiales-Clostridia
        Bacteria-Firmicutes-Clostridiales-Clostridia
        Bacteria-Firmicutes-Clostridiales-Clostridia
    /],
    'classifications match',
);

# compare metrics
my $metrics_file = $tmp_rdp_file.'.metrics';
my $fh = eval{ Genome::Sys->open_file_for_reading($tmp_rdp_file.'.metrics'); };
die "Failed to open metrics file ($metrics_file): $@" if not $fh;
my %metrics;
while ( my $line = $fh->getline ) {
    chomp $line;
    my ($key, $val) =split('=', $line);
    $metrics{$key} = $val;
}
$fh->close;
is_deeply(\%metrics, { total => 10, success => 10, error => 0, }, 'metrics match');

done_testing();
exit;

