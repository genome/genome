-- Deploy result_software_result_permission
-- requires: result_software_result

BEGIN;

REVOKE ALL ON TABLE result.software_result FROM PUBLIC;
REVOKE ALL ON TABLE result.software_result FROM genome;
GRANT ALL ON TABLE result.software_result TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE result.software_result TO "gms-user";

COMMIT;
