-- Verify timeline_allocation_event_type

BEGIN;

SELECT id
FROM timeline.allocation_event_type
WHERE FALSE;

ROLLBACK;
