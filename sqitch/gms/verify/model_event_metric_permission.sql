-- Verify model_event_metric_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.event_metric', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.event_metric', 'SELECT')::int;

ROLLBACK;
