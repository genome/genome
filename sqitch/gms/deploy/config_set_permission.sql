-- Deploy config_set_permission
-- requires: config_set

BEGIN;

REVOKE ALL ON TABLE config.set FROM PUBLIC;
REVOKE ALL ON TABLE config.set FROM genome;
GRANT ALL ON TABLE config.set TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.set TO "gms-user";

COMMIT;
