-- Revert config.instrument_data_analysis_project_bridge_index_fail_count

BEGIN;

DROP INDEX config.instrument_data_analysis_project_bridge_fail_count_idx;

COMMIT;
