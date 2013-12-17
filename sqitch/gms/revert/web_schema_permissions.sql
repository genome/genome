-- Revert web_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA web FROM "gms-user";

COMMIT;
