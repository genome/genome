-- Deploy model.model.subclass_name
-- requires: model_model

BEGIN;

CREATE INDEX model_subclass_index on model.model using btree (subclass_name);

COMMIT;
