package Genome::VariantReviewDetail;

use strict;
use warnings;

use Genome;
class Genome::VariantReviewDetail {
    type_name => 'variant review detail',
    table_name => 'VARIANT_REVIEW_DETAIL',
    id_by => [
        id => { is => 'NUMBER', len => 10 },
    ],
    has => [
        subject_name            => { is => 'VARCHAR2', len => 1024 },
        supporting_dbs          => { is => 'VARCHAR2', len => 1024, is_optional => 1 },
        supporting_samples      => { is => 'VARCHAR2', len => 1024, is_optional => 1 },
        genes                   => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        variant_type            => { is => 'VARCHAR2', len => 1, is_optional => 1 },
        chromosome              => { is => 'VARCHAR2', len => 32 },
        start_position          => { is => 'NUMBER', len => 10 },
        stop_position           => { is => 'NUMBER', len => 10 },
        variant_seq_id          => { is => 'NUMBER', len => 10, is_optional => 1 },
        roi_seq_id              => { is => 'NUMBER', len => 10, is_optional => 1 },
        rgg_id                  => { is => 'NUMBER', len => 10, is_optional => 1 },
        delete_sequence         => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        insert_sequence_allele1 => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        insert_sequence_allele2 => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        variant_length          => { is => 'NUMBER', len => 10, is_optional => 1 },
        sample_name             => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        project_type            => { is => 'VARCHAR2', len => 1, is_optional => 1 },
    ],
    has_many => [
        reviews => { is => 'Genome::SNVManualReview', reverse_id_by => 'detail_id' },
    ],
    unique_constraints => [
        { properties => [qw/chromosome sample_name start_position subject_name/], sql => 'VRD_SCS_U' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
