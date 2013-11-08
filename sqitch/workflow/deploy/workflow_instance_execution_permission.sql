-- Deploy workflow_instance_execution_permission
-- requires: workflow_instance_execution

BEGIN;

REVOKE ALL ON TABLE workflow.instance_execution FROM PUBLIC;
REVOKE ALL ON TABLE workflow.instance_execution FROM genome;
GRANT ALL ON TABLE workflow.instance_execution TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE workflow.instance_execution TO "gms-user";

COMMIT;
