-- Deploy config.instrument_data_analysis_project_bridge_index_instrument_data_id
-- requires: config_instrument_data_analysis_project_bridge

BEGIN;

CREATE INDEX c_idapb_index_instrument_data_id ON config.instrument_data_analysis_project_bridge (instrument_data_id);

COMMIT;
