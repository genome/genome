-- Revert workflow_execution_metric

BEGIN;

DROP TABLE workflow.execution_metric;

COMMIT;
