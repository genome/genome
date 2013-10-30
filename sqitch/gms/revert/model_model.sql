-- Revert model_model

BEGIN;

DROP TABLE IF EXISTS model.model;

COMMIT;
