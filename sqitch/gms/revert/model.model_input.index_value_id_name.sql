-- Revert model.model_input.index_value_id_name

BEGIN;

DROP INDEX model.model_input_value_id_name_index;

COMMIT;
