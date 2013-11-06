-- Deploy config_schema_permissions
-- requires: config_schema

BEGIN;

REVOKE ALL ON SCHEMA config FROM PUBLIC;
REVOKE ALL ON SCHEMA config FROM genome;
GRANT ALL ON SCHEMA config TO genome;
GRANT USAGE ON SCHEMA config TO "gms-user";

COMMIT;
