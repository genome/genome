-- Deploy model.processing_profile_param.name
-- requires: model_processing_profile_param

BEGIN;

CREATE INDEX processing_profile_param_name_index on model.processing_profile_param using btree (name);

COMMIT;
