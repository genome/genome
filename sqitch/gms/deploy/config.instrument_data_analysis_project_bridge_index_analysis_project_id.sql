-- Deploy config.instrument_data_analysis_project_bridge_index_analysis_project_id
-- requires: config_instrument_data_analysis_project_bridge

BEGIN;

CREATE INDEX c_idapb_index_analysis_project_id ON config.instrument_data_analysis_project_bridge (analysis_project_id);

COMMIT;
