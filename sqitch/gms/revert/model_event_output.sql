-- Revert model_event_output

BEGIN;

DROP TABLE IF EXISTS model.event_output;

COMMIT;
