-- Revert config.instrument_data_analysis_project_bridge_index_status

BEGIN;

DROP INDEX config.instrument_data_analysis_project_bridge_status_idx;

COMMIT;
