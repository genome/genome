-- Revert disk.allocation.index_allocation_path

BEGIN;

DROP INDEX disk.allocation_allocation_path_index;

COMMIT;
