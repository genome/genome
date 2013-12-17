-- Verify model_build_link_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.build_link', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.build_link', 'SELECT')::int;

ROLLBACK;
