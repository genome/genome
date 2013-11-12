-- Verify config_instrument_data_analysis_project_bridge

BEGIN;

SELECT id, instrument_data_id, analysis_project_id, created_at,
    updated_at, status, reason, fail_count
FROM config.instrument_data_analysis_project_bridge
WHERE FALSE;

ROLLBACK;
