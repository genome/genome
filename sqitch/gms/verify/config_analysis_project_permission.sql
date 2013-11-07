-- Verify config_analysis_project_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'config.analysis_project', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'config.analysis_project', 'SELECT')::int;

ROLLBACK;
