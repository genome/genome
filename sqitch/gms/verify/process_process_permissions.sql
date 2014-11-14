-- Verify process_process_permissions

BEGIN;

SELECT 1/(1 - has_table_privilege('gms-user', 'process.process', 'DELETE')::int);

ROLLBACK;
