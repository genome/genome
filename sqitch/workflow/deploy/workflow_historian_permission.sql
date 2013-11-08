-- Deploy workflow_historian_permission
-- requires: workflow_historian

BEGIN;

REVOKE ALL ON TABLE workflow.historian FROM PUBLIC;
REVOKE ALL ON TABLE workflow.historian FROM genome;
GRANT ALL ON TABLE workflow.historian TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE workflow.historian TO "gms-user";

COMMIT;
