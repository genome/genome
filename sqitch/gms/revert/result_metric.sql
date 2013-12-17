-- Revert result_metric

BEGIN;

DROP TABLE IF EXISTS result.metric;

COMMIT;
