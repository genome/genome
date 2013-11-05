-- Deploy model.processing_profile.name
-- requires: model_processing_profile

BEGIN;

CREATE INDEX processing_profile_name_index on model.processing_profile using btree (name);

COMMIT;
