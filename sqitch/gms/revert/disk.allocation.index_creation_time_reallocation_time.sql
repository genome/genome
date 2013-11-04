-- Revert disk.allocation.index_creation_time_reallocation_time

BEGIN;

DROP INDEX disk.allocation_creation_reallocation_time_index;

COMMIT;
