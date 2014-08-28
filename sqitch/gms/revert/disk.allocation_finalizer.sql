-- Revert disk.allocation_finalizer

BEGIN;

ALTER TABLE disk.allocation DROP COLUMN finalizer_id;
DROP TABLE disk.allocation_finalizer;

COMMIT;
