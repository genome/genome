-- Revert result_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA result FROM "gms-user";

COMMIT;
