-- Revert process_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA process FROM "gms-user";

COMMIT;
