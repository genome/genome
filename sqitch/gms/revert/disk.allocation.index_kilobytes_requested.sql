-- Revert disk.allocation.index_kilobytes_requested

BEGIN;

DROP INDEX disk.kilobytes_requested_index;

COMMIT;
