-- Revert model_event_input

BEGIN;

DROP TABLE IF EXISTS model.event_input;

COMMIT;
