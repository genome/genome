-- Deploy subject_role_permission
-- requires: subject_role

BEGIN;

REVOKE ALL ON TABLE subject.role FROM PUBLIC;
REVOKE ALL ON TABLE subject.role FROM genome;
GRANT ALL ON TABLE subject.role TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.role TO "gms-user";

COMMIT;
