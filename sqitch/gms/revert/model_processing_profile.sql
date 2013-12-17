-- Revert model_processing_profile

BEGIN;

DROP TABLE IF EXISTS model.processing_profile;

COMMIT;
