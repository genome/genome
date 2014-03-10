-- Verify model_event_output

BEGIN;

SELECT event_id, param_name, param_value
FROM model.event_output
WHERE FALSE;

ROLLBACK;
