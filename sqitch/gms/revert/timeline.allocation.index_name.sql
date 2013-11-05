-- Revert timeline.allocation.index_name

BEGIN;

DROP INDEX timeline.allocation_name_idx;

COMMIT;
