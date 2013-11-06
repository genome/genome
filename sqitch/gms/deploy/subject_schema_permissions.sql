-- Deploy subject_schema_permissions
-- requires: subject_schema

BEGIN;

REVOKE ALL ON SCHEMA subject FROM PUBLIC;
REVOKE ALL ON SCHEMA subject FROM genome;
GRANT ALL ON SCHEMA subject TO genome;
GRANT USAGE ON SCHEMA subject TO "gms-user";

COMMIT;
