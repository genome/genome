-- Verify disk_file_summary_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'disk.file_summary', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'disk.file_summary', 'SELECT')::int;

ROLLBACK;
