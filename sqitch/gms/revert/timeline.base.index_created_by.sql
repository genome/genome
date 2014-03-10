-- Revert timeline.base.index_created_by

BEGIN;

DROP INDEX timeline.base_created_by_idx;

COMMIT;
