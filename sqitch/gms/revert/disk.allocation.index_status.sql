-- Revert disk.allocation.index_status

BEGIN;

DROP INDEX disk.allocation_status_idx;

COMMIT;
