-- Deploy result_schema_permissions
-- requires: result_schema

BEGIN;

REVOKE ALL ON SCHEMA result FROM PUBLIC;
REVOKE ALL ON SCHEMA result FROM genome;
GRANT ALL ON SCHEMA result TO genome;
GRANT USAGE ON SCHEMA result TO "gms-user";

COMMIT;
