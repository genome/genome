-- Deploy result.input.input_value_input_name
-- requires: result_input

BEGIN;

CREATE INDEX sri_iniv on result.input using btree (input_value, input_name);

COMMIT;
