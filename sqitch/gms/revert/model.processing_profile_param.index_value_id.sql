-- Revert model.processing_profile_param.index_value_id

BEGIN;

DROP INDEX model.processing_profile_param_value_id_index;

COMMIT;
