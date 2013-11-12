-- Deploy config.instrument_data_analysis_project_bridge_index_fail_count
-- requires: config_instrument_data_analysis_project_bridge

BEGIN;

CREATE INDEX c_adapb_index_fail_count ON config.instrument_data_analysis_project_bridge (fail_count);

COMMIT;
