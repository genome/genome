-- Revert config.instrument_data_analysis_project_bridge_index_status

BEGIN;

DROP INDEX config.c_idapb_status_index;

COMMIT;
