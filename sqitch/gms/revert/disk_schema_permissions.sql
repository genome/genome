-- Revert disk_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA disk FROM "gms-user";

COMMIT;
