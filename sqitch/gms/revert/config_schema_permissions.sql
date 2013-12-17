-- Revert config_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA config FROM "gms-user";

COMMIT;
