-- Revert model.processing_profile_param.index_value_class_name_value_id

BEGIN;

DROP INDEX model.processing_profile_param_name_value_class_id_index;

COMMIT;
