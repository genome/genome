-- Deploy model.processing_profile_param.value_id
-- requires: model_processing_profile_param

BEGIN;

CREATE INDEX processing_profile_param_value_id_index on model.processing_profile_param using btree (value_id);

COMMIT;
