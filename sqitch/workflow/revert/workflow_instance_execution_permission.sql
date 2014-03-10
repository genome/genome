-- Revert workflow_instance_execution_permission

BEGIN;

REVOKE ALL ON TABLE workflow.instance_execution FROM "gms-user";

COMMIT;
