-- Deploy web_schema_permissions
-- requires: web_schema

BEGIN;

REVOKE ALL ON SCHEMA web FROM PUBLIC;
REVOKE ALL ON SCHEMA web FROM genome;
GRANT ALL ON SCHEMA web TO genome;
GRANT USAGE ON SCHEMA web TO "gms-user";

COMMIT;
