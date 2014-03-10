-- Revert workflow_plan_permission

BEGIN;

REVOKE ALL ON TABLE workflow.plan FROM "gms-user";

COMMIT;
