-- Deploy workflow_schema_permission
-- requires: workflow_schema

BEGIN;

REVOKE ALL ON SCHEMA workflow FROM PUBLIC;
REVOKE ALL ON SCHEMA workflow FROM genome;
GRANT ALL ON SCHEMA workflow TO genome;
GRANT USAGE ON SCHEMA workflow TO "gms-user";

COMMIT;
