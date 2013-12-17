-- Deploy config.instrument_data_analysis_project_bridge_index_status
-- requires: config_instrument_data_analysis_project_bridge

BEGIN;

CREATE INDEX instrument_data_analysis_project_bridge_status_idx ON config.instrument_data_analysis_project_bridge (status);

COMMIT;
