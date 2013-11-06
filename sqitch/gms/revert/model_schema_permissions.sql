-- Revert model_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA model FROM "gms-user";

COMMIT;
