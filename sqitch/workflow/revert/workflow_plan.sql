-- Revert workflow_plan

BEGIN;

DROP TABLE workflow.plan;

COMMIT;
