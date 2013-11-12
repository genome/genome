-- Revert config.instrument_data_analysis_project_bridge_index_analysis_project_id

BEGIN;

DROP INDEX config.c_idapb_index_analysis_project_id;

COMMIT;
