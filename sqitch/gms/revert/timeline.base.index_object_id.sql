-- Revert timeline.base.index_object_id

BEGIN;

DROP INDEX timeline.base_object_id_idx;

COMMIT;
