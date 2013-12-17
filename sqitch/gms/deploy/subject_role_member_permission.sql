-- Deploy subject_role_member_permission
-- requires: subject_role_member

BEGIN;

REVOKE ALL ON TABLE subject.role_member FROM PUBLIC;
REVOKE ALL ON TABLE subject.role_member FROM genome;
GRANT ALL ON TABLE subject.role_member TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject.role_member TO "gms-user";

COMMIT;
