-- Verify disk_volume_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'disk.volume', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'disk.volume', 'SELECT')::int;

ROLLBACK;
