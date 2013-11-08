-- Verify workflow_execution_metric

BEGIN;

SELECT workflow_execution_id, name, value
FROM workflow.execution_metric
WHERE FALSE;

ROLLBACK;
