-- Verify model_processing_profile_param

BEGIN;

SELECT processing_profile_id, param_name, param_value, name, value_class_name, value_id
FROM model.processing_profile_param
WHERE FALSE;

ROLLBACK;
