-- Verify model_build_metric_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.build_metric', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.build_metric', 'SELECT')::int;

ROLLBACK;
