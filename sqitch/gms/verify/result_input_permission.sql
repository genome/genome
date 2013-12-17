-- Verify result_input_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'result.input', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'result.input', 'SELECT')::int;

ROLLBACK;
