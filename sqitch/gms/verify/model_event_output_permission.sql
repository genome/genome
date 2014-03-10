-- Verify model_event_output_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.event_output', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.event_output', 'SELECT')::int;

ROLLBACK;
