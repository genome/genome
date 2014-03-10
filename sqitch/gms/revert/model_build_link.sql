-- Revert model_build_link

BEGIN;

DROP TABLE IF EXISTS model.build_link;

COMMIT;
