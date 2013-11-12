-- Revert config.instrument_data_analysis_project_bridge_index_fail_count

BEGIN;

DROP INDEX config.c_adapb_index_fail_count;

COMMIT;
