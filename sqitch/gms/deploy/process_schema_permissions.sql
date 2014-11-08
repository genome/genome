-- Deploy process_schema_permissions
-- requires: process_schema

BEGIN;

REVOKE ALL ON SCHEMA process FROM PUBLIC;
REVOKE ALL ON SCHEMA process FROM genome;
GRANT ALL ON SCHEMA process TO genome;
GRANT USAGE ON SCHEMA process TO "gms-user";

COMMIT;
