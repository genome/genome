-- Revert model.model.created_by_idx

BEGIN;

DROP INDEX model.model__created_by_idx;

COMMIT;
