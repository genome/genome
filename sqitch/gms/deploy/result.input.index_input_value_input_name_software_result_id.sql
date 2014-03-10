-- Deploy result.input.input_value_input_name_software_result_id
-- requires: result_input

BEGIN;

CREATE INDEX sri_inivsri on result.input using btree (input_value, input_name, software_result_id);

COMMIT;
