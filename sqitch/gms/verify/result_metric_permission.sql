-- Verify result_metric_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'result.metric', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'result.metric', 'SELECT')::int;

ROLLBACK;
