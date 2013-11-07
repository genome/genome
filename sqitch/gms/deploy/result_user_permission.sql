-- Deploy result_user_permission
-- requires: result_user

BEGIN;

REVOKE ALL ON TABLE result."user" FROM PUBLIC;
REVOKE ALL ON TABLE result."user" FROM genome;
GRANT ALL ON TABLE result."user" TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE result."user" TO "gms-user";

COMMIT;
