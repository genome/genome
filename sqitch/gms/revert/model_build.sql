-- Revert model_build

BEGIN;

DROP TABLE IF EXISTS model.build;

COMMIT;
