-- Deploy subject_user_permission
-- requires: subject_user

BEGIN;

REVOKE ALL ON TABLE subject."user" FROM PUBLIC;
REVOKE ALL ON TABLE subject."user" FROM genome;
GRANT ALL ON TABLE subject."user" TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE subject."user" TO "gms-user";

COMMIT;
