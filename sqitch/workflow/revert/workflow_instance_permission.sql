-- Revert workflow_instance_permission

BEGIN;

REVOKE ALL ON TABLE workflow.instance FROM "gms-user";

COMMIT;
