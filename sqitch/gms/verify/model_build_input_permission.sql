-- Verify model_build_input_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.build_input', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.build_input', 'SELECT')::int;

ROLLBACK;
