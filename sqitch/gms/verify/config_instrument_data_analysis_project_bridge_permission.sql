-- Verify config_instrument_data_analysis_project_bridge_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'config.instrument_data_analysis_project_bridge', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'config.instrument_data_analysis_project_bridge', 'SELECT')::int;

ROLLBACK;
