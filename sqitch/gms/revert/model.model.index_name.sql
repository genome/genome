-- Revert model.model.index_name

BEGIN;

DROP INDEX model.model_name_index;

COMMIT;
