-- Verify model_build_metric

BEGIN;

SELECT build_id, metric_name, metric_value
FROM model.build_metric
WHERE FALSE;

ROLLBACK;
