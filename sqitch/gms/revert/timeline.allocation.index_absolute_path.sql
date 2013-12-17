-- Revert timeline.allocation.index_absolute_path

BEGIN;

DROP INDEX timeline.allocation_absolute_path_idx;

COMMIT;
