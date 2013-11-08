-- Deploy workflow_service_permission
-- requires: workflow_service

BEGIN;

REVOKE ALL ON TABLE workflow.service FROM PUBLIC;
REVOKE ALL ON TABLE workflow.service FROM genome;
GRANT ALL ON TABLE workflow.service TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE workflow.service TO "gms-user";

COMMIT;
