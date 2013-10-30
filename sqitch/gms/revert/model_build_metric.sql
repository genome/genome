-- Revert model_build_metric

BEGIN;

DROP TABLE IF EXISTS model.build_metric;

COMMIT;
