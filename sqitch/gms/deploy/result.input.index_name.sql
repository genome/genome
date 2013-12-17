-- Deploy result.input.name
-- requires: result_input

BEGIN;

CREATE INDEX result_input_name_index on result.input using btree (name);

COMMIT;
