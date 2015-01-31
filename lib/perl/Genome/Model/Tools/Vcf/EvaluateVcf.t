#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Vcf::EvaluateVcf;

my $vcf = "/gscmnt/gc2801/analytics/tabbott/vcf-evaluate/src/vcf-evaluation/tmp/small.vcf.gz";
#/gscmnt/gc13028/info/model_data/6195f53ff81046959a3bfee124c186f6/buildf60b23fdb26c45f486cf2d5937a7502b/variants/snvs.vcf.gz";

my $cmd = Genome::Model::Tools::Vcf::EvaluateVcf->create(
        vcf => $vcf,
        old_sample => "H_IJ-NA12878",
        new_sample => "NA12878",
        gold_vcf => "/gscmnt/gc2801/analytics/tabbott/vcf-evaluate/gold/snvs.vcf.gz",
        gold_sample => "NA12878",
        roi => "/gscmnt/gc2801/analytics/tabbott/vcf-evaluate/gold/highconf.bed.gz",
        true_negative_bed => "/gscmnt/gc2801/analytics/tabbott/vcf-evaluate/gold/tn.bed.gz",
        output_directory => "/gscmnt/gc2801/analytics/idas/vcf/test",
    );

$cmd->execute;
