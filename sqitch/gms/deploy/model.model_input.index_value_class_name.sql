-- Deploy model.model_input.value_class_name
-- requires: model_model_input

BEGIN;

CREATE INDEX model_input_value_class_index on model.model_input using btree (value_class_name);

COMMIT;
