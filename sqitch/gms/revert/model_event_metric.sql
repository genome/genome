-- Revert model_event_metric

BEGIN;

DROP TABLE IF EXISTS model.event_metric;

COMMIT;
