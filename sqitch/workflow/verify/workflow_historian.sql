-- Verify workflow_historian

BEGIN;

SELECT net_key, operation_id, color, workflow_instance_id
FROM workflow.historian
WHERE FALSE;

ROLLBACK;
