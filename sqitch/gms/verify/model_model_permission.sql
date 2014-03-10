-- Verify model_model_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.model', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.model', 'SELECT')::int;

ROLLBACK;
