-- Revert result.input.index_software_result_id_name

BEGIN;

DROP INDEX result.result_input_id_name;

COMMIT;
