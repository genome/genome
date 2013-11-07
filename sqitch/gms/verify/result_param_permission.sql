-- Verify result_param_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'result.param', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'result.param', 'SELECT')::int;

ROLLBACK;
