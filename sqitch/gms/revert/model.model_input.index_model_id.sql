-- Revert model.model_input.index_model_id

BEGIN;

DROP INDEX model.model_input_model_id_index;

COMMIT;
