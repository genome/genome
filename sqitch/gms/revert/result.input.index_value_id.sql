-- Revert result.input.index_value_id

BEGIN;

DROP INDEX result.result_input_value_id_index;

COMMIT;
