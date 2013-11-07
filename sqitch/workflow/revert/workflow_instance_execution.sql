-- Revert workflow_instance_execution

BEGIN;

DROP TABLE workflow.instance_execution;

COMMIT;
