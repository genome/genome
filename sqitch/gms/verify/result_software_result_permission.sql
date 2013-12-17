-- Verify result_software_result_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'result.software_result', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'result.software_result', 'SELECT')::int;

ROLLBACK;
