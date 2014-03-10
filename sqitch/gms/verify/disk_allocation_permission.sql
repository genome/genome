-- Verify disk_allocation_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'disk.allocation', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'disk.allocation', 'SELECT')::int;

ROLLBACK;
