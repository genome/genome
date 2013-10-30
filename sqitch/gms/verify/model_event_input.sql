-- Verify model_event_input

BEGIN;

SELECT event_id, param_name, param_value
FROM model.event_input
WHERE FALSE;

ROLLBACK;
