package Genome::Capture::Target;

use strict;
use warnings;

use Genome;

class Genome::Capture::Target {
    table_name => q|
        (select
            AMPLIFICATION_TARGET_ID id,
            AMPLIFICATION_ROI_ID region_id,
            TARGET_STAG_ID tag_id,
            CREATION_EVENT_ID pse_id
        from amplification_target
        ) capture_target
    |,
    id_by => [
        id => { },
    ],
    has => [
        region_id => { },
        region => {
            is => 'Genome::Capture::Region',
            id_by => 'region_id',
        },
        tag_id => { },
        pse_id => { },
    ],
    doc         => 'the target, usually exon or cds',
    data_source => 'Genome::DataSource::Oltp',
};

1;
