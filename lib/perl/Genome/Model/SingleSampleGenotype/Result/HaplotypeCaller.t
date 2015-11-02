#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::InstrumentData::AlignmentResult::Merged::Speedseq;
use Genome::Test::Factory::Model::ImportedReferenceSequence;

use File::Spec;
use Sub::Override;

use Test::More tests => 2;

my $pkg = 'Genome::Model::SingleSampleGenotype::Result::HaplotypeCaller';
use_ok($pkg);

my $reference_model = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object;
my $reference = Genome::Test::Factory::Build->setup_object(model_id => $reference_model->id);
Genome::Sys->write_file(File::Spec->join($reference->data_directory, 'all_sequences.fa'), '>null');

my $speedseq_result = Genome::Test::Factory::InstrumentData::AlignmentResult::Merged::Speedseq->setup_object(
#    instrument_data => [$build->instrument_data],
    reference_build => $reference,
#    map {; $_ => $model->$_ } (qw(aligner_name aligner_version aligner_params)),
);

require Genome::Model::Tools::Gatk::HaplotypeCaller;
my $override = Sub::Override->new(
    'Genome::Model::Tools::Gatk::Base::_execute_body',
    sub {
        my $self = shift;
        Genome::Sys->write_file($self->output_vcf, '## NO DATA');
        return 1;
    }
);

my $result = $pkg->create(
    alignment_result => $speedseq_result,
    intervals => [1],
    haplotype_caller_version => '3.4',
    emit_reference_confidence => 'GVCF',
);
isa_ok($result, $pkg, 'generated result');

