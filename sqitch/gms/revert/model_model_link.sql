-- Revert model_model_link

BEGIN;

DROP TABLE IF EXISTS model.model_link;

COMMIT;
