-- Deploy subject_project_permission
-- requires: subject_project

BEGIN;

REVOKE ALL ON TABLE subject.project FROM PUBLIC;
REVOKE ALL ON TABLE subject.project FROM genome;
GRANT ALL ON TABLE subject.project TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.project TO "gms-user";

COMMIT;
