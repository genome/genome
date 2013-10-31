-- Revert timeline_allocation

BEGIN;

DROP TABLE IF EXISTS timeline.allocation;

COMMIT;
