-- Revert workflow_service_permission

BEGIN;

REVOKE ALL ON TABLE workflow.service FROM "gms-user";

COMMIT;
