-- Revert model.model.index_subclass_name

BEGIN;

DROP INDEX model.model_subclass_index;

COMMIT;
