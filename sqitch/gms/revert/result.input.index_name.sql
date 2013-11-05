-- Revert result.input.index_name

BEGIN;

DROP INDEX result.result_input_name_index;

COMMIT;
