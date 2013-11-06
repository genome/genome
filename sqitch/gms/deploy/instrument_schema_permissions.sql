-- Deploy instrument_schema_permissions
-- requires: instrument_schema

BEGIN;

REVOKE ALL ON SCHEMA instrument FROM PUBLIC;
REVOKE ALL ON SCHEMA instrument FROM genome;
GRANT ALL ON SCHEMA instrument TO genome;
GRANT USAGE ON SCHEMA instrument TO "gms-user";

COMMIT;
