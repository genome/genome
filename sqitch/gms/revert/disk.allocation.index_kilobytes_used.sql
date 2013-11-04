-- Revert disk.allocation.index_kilobytes_used

BEGIN;

DROP INDEX disk.kilobytes_used_index;

COMMIT;
