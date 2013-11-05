-- Revert model.processing_profile_param.index_name

BEGIN;

DROP INDEX model.processing_profile_param_name_index;

COMMIT;
