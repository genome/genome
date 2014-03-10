-- Deploy model_schema_permissions
-- requires: model_schema

BEGIN;

REVOKE ALL ON SCHEMA model FROM PUBLIC;
REVOKE ALL ON SCHEMA model FROM genome;
GRANT ALL ON SCHEMA model TO genome;
GRANT USAGE ON SCHEMA model TO "gms-user";

COMMIT;
