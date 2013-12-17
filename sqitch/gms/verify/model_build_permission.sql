-- Verify model_build_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.build', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.build', 'SELECT')::int;

ROLLBACK;
