-- Deploy model.model_input.value_id_name
-- requires: model_model_input

BEGIN;

CREATE INDEX model_input_value_id_name_index on model.model_input using btree (value_id, name);

COMMIT;
