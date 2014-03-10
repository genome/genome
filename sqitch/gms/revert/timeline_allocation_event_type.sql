-- Revert timeline_allocation_event_type

BEGIN;

DROP TABLE IF EXISTS timeline.allocation_event_type;

COMMIT;
