-- Revert model_event

BEGIN;

DROP TABLE IF EXISTS model.event;

COMMIT;
