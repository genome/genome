-- Deploy timeline_schema_permissions
-- requires: timeline_schema

BEGIN;

REVOKE ALL ON SCHEMA timeline FROM PUBLIC;
REVOKE ALL ON SCHEMA timeline FROM genome;
GRANT ALL ON SCHEMA timeline TO genome;
GRANT USAGE ON SCHEMA timeline TO "gms-user";

COMMIT;
