-- Deploy result.input.value_id
-- requires: result_input

BEGIN;

CREATE INDEX result_input_value_id_index on result.input using btree (value_id);

COMMIT;
