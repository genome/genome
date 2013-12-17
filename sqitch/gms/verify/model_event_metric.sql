-- Verify model_event_metric

BEGIN;

SELECT event_id, metric_name, metric_value
FROM model.event_metric
WHERE FALSE;

ROLLBACK;
