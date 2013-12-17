-- Revert result.input.index_input_value_input_name_software_result_id

BEGIN;

DROP INDEX result.sri_inivsri;

COMMIT;
