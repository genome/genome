-- Revert model_model_group

BEGIN;

DROP TABLE IF EXISTS model.model_group;

COMMIT;
