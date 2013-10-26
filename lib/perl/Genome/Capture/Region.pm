package Genome::Capture::Region;

use strict;
use warnings;

use Genome;

class Genome::Capture::Region {
    table_name => q|
        (select
            AMPLIFICATION_ROI_ID id,
            HTMP_PROJECT_ID project_id,
            REGION_OF_INTEREST_NAME name,
            REGION_OF_INTEREST_SEQ_ID seq_id,
            CREATION_EVENT_ID pse_id
        from amplification_roi
        ) capture_region
    |,
    id_by => [
        id => { },
    ],
    has => [
        project_id => { },
        name => { },
        seq_id => { },
        pse_id => { },
    ],
    doc         => 'usually a gene or some other region of interest',
    data_source => 'Genome::DataSource::Oltp',
};


1;
