-- Verify disk_group_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'disk."group"', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'disk."group"', 'SELECT')::int;

ROLLBACK;
