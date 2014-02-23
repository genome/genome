-- Revert model.model.run_as_idx

BEGIN;

DROP INDEX model.model__run_as_idx;

COMMIT;
