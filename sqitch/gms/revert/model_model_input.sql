-- Revert model_model_input

BEGIN;

DROP TABLE IF EXISTS model.model_input;

COMMIT;
