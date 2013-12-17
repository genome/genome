-- Revert model.model.index_processing_profile_id

BEGIN;

DROP INDEX model.model_processing_profile_index;

COMMIT;
