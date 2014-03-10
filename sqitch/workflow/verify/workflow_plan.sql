-- Verify workflow_plan

BEGIN;

SELECT workflow_plan_id, xml
FROM workflow.plan
WHERE FALSE;

ROLLBACK;
