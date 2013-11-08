-- Deploy workflow_instance_permission
-- requires: workflow_instance

BEGIN;

REVOKE ALL ON TABLE workflow.instance FROM PUBLIC;
REVOKE ALL ON TABLE workflow.instance FROM genome;
GRANT ALL ON TABLE workflow.instance TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE workflow.instance TO "gms-user";

COMMIT;
