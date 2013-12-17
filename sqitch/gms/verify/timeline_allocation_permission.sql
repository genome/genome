-- Verify timeline_allocation_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'timeline.allocation', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'timeline.allocation', 'SELECT')::int;

ROLLBACK;
