-- Revert timeline_base

BEGIN;

DROP TABLE IF EXISTS timeline.base;

COMMIT;
