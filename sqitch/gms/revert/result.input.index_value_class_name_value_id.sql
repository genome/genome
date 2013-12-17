-- Revert result.input.index_value_class_name_value_id

BEGIN;

DROP INDEX result.result_input_value_class_id_index;

COMMIT;
