-- Revert workflow_historian

BEGIN;

DROP TABLE workflow.historian;

COMMIT;
