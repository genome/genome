-- Verify timeline_base_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'timeline.base', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'timeline.base', 'SELECT')::int;

ROLLBACK;
