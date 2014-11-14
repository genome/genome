-- Verify process_status_event_permissions

BEGIN;

SELECT 1/(1 - has_table_privilege('gms-user', 'process.status_event', 'DELETE')::int);
SELECT 1/(1 - has_table_privilege('gms-user', 'process.status_event', 'UPDATE')::int);

ROLLBACK;
