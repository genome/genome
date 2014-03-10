-- Revert model.model_input.index_value_class_name

BEGIN;

DROP INDEX model.model_input_value_class_index;

COMMIT;
