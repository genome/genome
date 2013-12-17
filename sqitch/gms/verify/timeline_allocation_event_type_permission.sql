-- Verify timeline_allocation_event_type_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'timeline.allocation_event_type', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'timeline.allocation_event_type', 'SELECT')::int;

ROLLBACK;
