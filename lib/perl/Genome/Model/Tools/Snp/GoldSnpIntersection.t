#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 5;
use File::Slurp;

use above 'Genome';

BEGIN {
        use_ok('Genome::Model::Tools::Snp::GoldSnpIntersection');
    };

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Snp/GoldSnpIntersection';

my $gold_snp_file = "$dir/gold.snp";
my $maq_snp_file  = "$dir/maq.snp";
my $sam_snp_file  = "$dir/sam.snp";

my $exp_maq_report = read_file($dir.'/report.maq.ori');
my $exp_sam_report = read_file($dir.'/report.sam.ori');

=cut

my $maq_gsi = Genome::Model::Tools::Snp::GoldSnpIntersection->create(
    snp_file      => $maq_snp_file,
    gold_snp_file => $gold_snp_file,
);

isa_ok($maq_gsi,'Genome::Model::Tools::Snp::GoldSnpIntersection');

=cut

my $maq_cmd = Genome::Model::Tools::Snp::GoldSnpIntersection->create(
    snp_file => $maq_snp_file,
    gold_snp_file => $gold_snp_file,
);
$maq_cmd->execute;
ok($maq_cmd->_report_txt, 'MAQ gold-snp-intersection execute');
is($maq_cmd->_report_txt, $exp_maq_report, 'MAQ gold-snp-intersection output matches the expected original one.')
    #or do {
    #        IO::File->new(">t.new")->print($maq_cmd->_report_txt);
    #        IO::File->new(">t.old")->print($exp_maq_report);
    #        system "kdiff3 t.old t.new";
    #    }
    ;

=cut

my $sam_gsi = Genome::Model::Tools::Snp::GoldSnpIntersection->create(
    snp_file      => $sam_snp_file,
    gold_snp_file => $gold_snp_file,
    sam_format    => 1,
);

my $sam_report = $sam_gsi->execute;

=cut

my $sam_cmd = Genome::Model::Tools::Snp::GoldSnpIntersection->create(
    snp_file => $sam_snp_file,
    gold_snp_file => $gold_snp_file,
    snp_format => "sam",
);
$sam_cmd->execute;
ok($sam_cmd->_report_txt,'SAM gold-snp-intersection execute ok');
is($sam_cmd->_report_txt, $exp_sam_report, 'SAM gold-snp-intersection output matches the expected original one.');
