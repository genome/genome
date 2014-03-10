-- Verify model_feature_list_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.feature_list', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.feature_list', 'SELECT')::int;

ROLLBACK;
