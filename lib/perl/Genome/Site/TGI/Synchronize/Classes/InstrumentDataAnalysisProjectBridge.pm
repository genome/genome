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
        swo.analysis_project_id
    FROM
    (
        SELECT to_char(i.analysis_id) instrument_data_id FROM index_illumina i
        UNION ALL
        SELECT to_char(ri454.seq_id) instrument_data_id FROM region_index_454 ri454
        UNION ALL
        SELECT to_char(g.seq_id) instrument_data_id FROM external_genotyping g
        UNION ALL
        SELECT to_char(g.seq_id) instrument_data_id FROM illumina_genotyping g
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

sub entity_name { return 'analysis project instrument data'; }

sub genome_class_for_create { return 'Genome::Config::AnalysisProject::InstrumentDataBridge'; }

sub genome_class_for_comparison { return 'Genome::Site::TGI::Synchronize::Classes::AnalysisProjectInstrumentData'; }

sub properties_to_copy {
    return ('instrument_data_id', 'analysis_project_id');
}

1;

