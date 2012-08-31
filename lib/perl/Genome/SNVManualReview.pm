package Genome::SNVManualReview;

use strict;
use warnings;

use Genome;
class Genome::SNVManualReview {
    type_name => 'snv manual review',
    table_name => 'SNV_MANUAL_REVIEW',
    id_by => [
        id => { is => 'NUMBER', len => 10 },
    ],
    has => [
        review_type                 => { is => 'VARCHAR2', len => 1, is_optional => 1 },
        reviewer                    => { is => 'VARCHAR2', len => 256 },
        notes                       => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        detail_id                   => { is => 'NUMBER', len => 10 },
        genotype_iub_code           => { is => 'VARCHAR2', len => 1 },
        pass_manual_review          => { is => 'VARCHAR2', len => 1 },
        manual_genotype_iub_normal  => { is => 'VARCHAR2', len => 1, is_optional => 1 },
        manual_genotype_iub_tumor   => { is => 'VARCHAR2', len => 1, is_optional => 1 },
        manual_genotype_iub_relapse => { is => 'VARCHAR2', len => 1, is_optional => 1 },
        somatic_status              => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        data_needed                 => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        build_id                    => { is => 'NUMBER', len => 10 },
        dump_date                   => { is => 'DATE', len => 19 },
    ],
    unique_constraints => [
        { properties => [qw/build_id detail_id dump_date reviewer/], sql => 'SMR_DBDR_U' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
