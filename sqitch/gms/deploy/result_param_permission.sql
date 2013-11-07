-- Deploy result_param_permission
-- requires: result_param

BEGIN;

REVOKE ALL ON TABLE result.param FROM PUBLIC;
REVOKE ALL ON TABLE result.param FROM genome;
GRANT ALL ON TABLE result.param TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE result.param TO "gms-user";

COMMIT;
