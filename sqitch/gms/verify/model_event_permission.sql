-- Verify model_event_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.event', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.event', 'SELECT')::int;

ROLLBACK;
