-- Revert model.processing_profile_param.index_processing_profile_id_param_name

BEGIN;

DROP INDEX model.processing_profile_param_id_param_name_index;

COMMIT;
