-- Revert workflow_instance

BEGIN;

DROP TABLE workflow.instance;

COMMIT;
