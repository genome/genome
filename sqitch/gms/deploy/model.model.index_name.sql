-- Deploy model.model.name
-- requires: model_model

BEGIN;

CREATE INDEX model_name_index on model.model using btree (name);

COMMIT;
