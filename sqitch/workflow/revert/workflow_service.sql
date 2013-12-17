-- Revert workflow_service

BEGIN;

DROP TABLE workflow.service;

COMMIT;
