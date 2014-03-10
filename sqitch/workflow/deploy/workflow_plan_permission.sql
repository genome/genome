-- Deploy workflow_plan_permission
-- requires: workflow_plan

BEGIN;

REVOKE ALL ON TABLE workflow.plan FROM PUBLIC;
REVOKE ALL ON TABLE workflow.plan FROM genome;
GRANT ALL ON TABLE workflow.plan TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE workflow.plan TO "gms-user";

COMMIT;
