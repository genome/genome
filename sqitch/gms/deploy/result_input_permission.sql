-- Deploy result_input_permission
-- requires: result_input

BEGIN;

REVOKE ALL ON TABLE result.input FROM PUBLIC;
REVOKE ALL ON TABLE result.input FROM genome;
GRANT ALL ON TABLE result.input TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE result.input TO "gms-user";

COMMIT;
