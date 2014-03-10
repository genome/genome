-- Revert timeline_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA timeline FROM "gms-user";

COMMIT;
