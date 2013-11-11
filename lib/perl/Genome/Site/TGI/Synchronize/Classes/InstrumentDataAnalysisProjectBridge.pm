package Genome::Site::TGI::Synchronize::Classes::InstrumentDataAnalysisProjectBridge;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::InstrumentDataAnalysisProjectBridge {
    table_name => <<'EOS'
(
    SELECT
        x.instrument_data_id instrument_data_id,
        swo.analysis_project_id
    FROM
    (
        SELECT to_char(i.analysis_id) instrument_data_id FROM index_illumina i
        UNION ALL
        SELECT to_char(g.seq_id) instrument_data_id FROM external_genotyping g
        UNION ALL
        SELECT to_char(g.seq_id) instrument_data_id FROM illumina_genotyping g WHERE g.status = 'PASS'
    ) x
    LEFT JOIN GSC.woi_sequence_product wsp ON wsp.seq_id = x.instrument_data_id
    LEFT JOIN work_order_item@oltp woi ON wsp.woi_id = woi.woi_id
    LEFT JOIN setup_work_order@oltp swo ON swo.setup_wo_id = woi.setup_wo_id
    WHERE swo.analysis_project_id IS NOT NULL
)
idapp
EOS
    ,
    id_by => [
        instrument_data_id => {
            is => 'Text',
        },
        analysis_project_id => {
            is => 'Text',
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub properties_to_copy {
    return ('instrument_data_id', 'analysis_project_id');
}

sub sync_id {
    my $self = shift;
    return join("\t", $self->instrument_data_id, $self->analysis_project_id);
}

1;
