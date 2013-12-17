-- Deploy model.processing_profile_param.value_class_name_value_id
-- requires: model_processing_profile_param

BEGIN;

CREATE INDEX processing_profile_param_name_value_class_id_index on model.processing_profile_param using btree (value_class_name, value_id);

COMMIT;
