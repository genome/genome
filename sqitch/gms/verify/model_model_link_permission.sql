-- Verify model_model_link_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.model_link', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.model_link', 'SELECT')::int;

ROLLBACK;
