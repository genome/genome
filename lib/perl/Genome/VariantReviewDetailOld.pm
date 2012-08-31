# review gsanders
# accordding to charris and dlarson this table was a hack for dlarson to see old data when jschindl was putting in new variant review stuff
# this can me removed and the table dropped IF AND ONLY IF we can show that the new variant review tables have inherited the information in this old table

package Genome::VariantReviewDetailOld;

use strict;
use warnings;

use Genome;
class Genome::VariantReviewDetailOld {
    type_name => 'variant review detail old',
    table_name => 'VARIANT_REVIEW_DETAIL_OLD',
    id_by => [
        id => { is => 'NUMBER', len => 10 },
    ],
    has => [
        begin_position          => { is => 'NUMBER', len => 10, is_optional => 1 },
        chromosome              => { is => 'VARCHAR2', len => 32, is_optional => 1 },
        data_needed             => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        delete_sequence         => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        end_position            => { is => 'NUMBER', len => 10, is_optional => 1 },
        finisher_3730_review    => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        finisher_manual_review  => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        genes                   => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        insert_sequence_allele1 => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        insert_sequence_allele2 => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        manual_genotype_normal  => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        manual_genotype_relapse => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        manual_genotype_tumor   => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        notes                   => { is => 'VARCHAR2', len => 4000, is_optional => 1 },
        pass_manual_review      => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        rgg_id                  => { is => 'NUMBER', len => 10, is_optional => 1 },
        roi_seq_id              => { is => 'NUMBER', len => 10, is_optional => 1 },
        sample_name             => { is => 'VARCHAR2', len => 1024, is_optional => 1 },
        somatic_status          => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        supporting_dbs          => { is => 'VARCHAR2', len => 1024, is_optional => 1 },
        supporting_samples      => { is => 'VARCHAR2', len => 1024, is_optional => 1 },
        variant_length          => { is => 'NUMBER', len => 10, is_optional => 1 },
        variant_seq_id          => { is => 'NUMBER', len => 10, is_optional => 1 },
        variant_type            => { is => 'VARCHAR2', len => 1, is_optional => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
