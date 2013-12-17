-- Deploy model.model_input.model_id
-- requires: model_model_input

BEGIN;

CREATE INDEX model_input_model_id_index on model.model_input using btree (model_id);

COMMIT;
