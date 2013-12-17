-- Deploy result.input.software_result_id_name
-- requires: result_input

BEGIN;

CREATE INDEX result_input_id_name on result.input using btree (software_result_id, name);

COMMIT;
