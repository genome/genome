-- Verify disk_volume_group_bridge_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'disk.volume_group_bridge', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'disk.volume_group_bridge', 'SELECT')::int;

ROLLBACK;
