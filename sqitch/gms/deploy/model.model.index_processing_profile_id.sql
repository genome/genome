-- Deploy model.model.processing_profile_id
-- requires: model_model

BEGIN;

CREATE INDEX model_processing_profile_index on model.model using btree (processing_profile_id);

COMMIT;
