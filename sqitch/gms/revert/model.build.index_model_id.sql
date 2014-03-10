-- Revert model.build.index_model_id

BEGIN;

DROP INDEX model.build_model_index;

COMMIT;
