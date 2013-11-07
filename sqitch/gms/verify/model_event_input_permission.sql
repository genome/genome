-- Verify model_event_input_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'model.event_input', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'model.event_input', 'SELECT')::int;

ROLLBACK;
