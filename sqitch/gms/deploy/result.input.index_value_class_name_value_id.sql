-- Deploy result.input.value_class_name_value_id
-- requires: result_input

BEGIN;

CREATE INDEX result_input_value_class_id_index on result.input using btree (value_class_name, value_id);

COMMIT;
