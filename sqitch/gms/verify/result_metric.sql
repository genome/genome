-- Verify result_metric

BEGIN;

SELECT software_result_id, metric_name, metric_value
FROM result.metric
WHERE FALSE;

ROLLBACK;
