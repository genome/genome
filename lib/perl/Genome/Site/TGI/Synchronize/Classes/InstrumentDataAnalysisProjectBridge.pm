package Genome::Site::TGI::Synchronize::Classes::InstrumentDataAnalysisProjectBridge;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::InstrumentDataAnalysisProjectBridge {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => <<'EOS'
(
    SELECT
        x.instrument_data_id instrument_data_id,
        x.analysis_project_id analysis_project_id
    FROM
    (
        SELECT i.analysis_id instrument_data_id, swo.analysis_project_id FROM index_illumina i
        LEFT JOIN GSC.woi_sequence_product wsp ON wsp.seq_id = i.seq_id
        LEFT JOIN work_order_item woi ON wsp.woi_id = woi.woi_id
        LEFT JOIN setup_work_order swo ON swo.setup_wo_id = woi.setup_wo_id
        WHERE swo.analysis_project_id IS NOT NULL
        UNION ALL
        SELECT g.seq_id instrument_data_id, swo.analysis_project_id FROM external_genotyping g
        LEFT JOIN GSC.woi_sequence_product wsp ON wsp.seq_id = g.seq_id
        LEFT JOIN work_order_item woi ON wsp.woi_id = woi.woi_id
        LEFT JOIN setup_work_order swo ON swo.setup_wo_id = woi.setup_wo_id
        WHERE swo.analysis_project_id IS NOT NULL
        UNION ALL
        SELECT g.seq_id instrument_data_id, swo.analysis_project_id FROM illumina_genotyping g
        LEFT JOIN GSC.woi_sequence_product wsp ON wsp.seq_id = g.seq_id
        LEFT JOIN work_order_item woi ON wsp.woi_id = woi.woi_id
        LEFT JOIN setup_work_order swo ON swo.setup_wo_id = woi.setup_wo_id
        WHERE swo.analysis_project_id IS NOT NULL
    ) x
)
idapp
EOS
    ,
    id_by => [
        instrument_data_id => {
            is => 'Number',
        },
        analysis_project_id => {
            is => 'Number',
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub entity_name { return 'analysis project instrument data'; }

sub genome_class_for_create { return 'Genome::Config::AnalysisProject::InstrumentDataBridge'; }

sub genome_class_for_comparison { return 'Genome::Site::TGI::Synchronize::Classes::AnalysisProjectInstrumentData'; }

sub properties_to_copy {
    return ('instrument_data_id', 'analysis_project_id');
}

1;

