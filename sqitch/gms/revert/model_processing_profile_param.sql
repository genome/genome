-- Revert model_processing_profile_param

BEGIN;

DROP TABLE IF EXISTS model.processing_profile_param;

COMMIT;
