-- Revert workflow_historian_permission

BEGIN;

REVOKE ALL ON TABLE workflow.historian FROM "gms-user";

COMMIT;
