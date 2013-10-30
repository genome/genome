-- Revert model_build_input

BEGIN;

DROP TABLE IF EXISTS model.build_input;

COMMIT;
