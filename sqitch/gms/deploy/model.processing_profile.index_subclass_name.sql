-- Deploy model.processing_profile.subclass_name
-- requires: model_processing_profile

BEGIN;

CREATE INDEX processing_profile_subclass_index on model.processing_profile using btree (subclass_name);

COMMIT;
