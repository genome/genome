-- Deploy model.processing_profile_param.processing_profile_id_name
-- requires: model_processing_profile_param

BEGIN;

CREATE INDEX processing_profile_param_id_name_index on model.processing_profile_param using btree (processing_profile_id, name);

COMMIT;
