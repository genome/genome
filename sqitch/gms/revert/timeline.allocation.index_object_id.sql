-- Revert timeline.allocation.index_object_id

BEGIN;

DROP INDEX timeline.allocation_object_id_idx;

COMMIT;
