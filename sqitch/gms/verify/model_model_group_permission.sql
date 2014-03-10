-- Verify model_model_group_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.model_group', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.model_group', 'SELECT')::int;

ROLLBACK;
