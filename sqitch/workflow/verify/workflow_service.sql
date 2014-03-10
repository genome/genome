-- Verify workflow_service

BEGIN;

SELECT hostname, port, start_time, process_id, username
FROM workflow.service
WHERE FALSE;

ROLLBACK;
