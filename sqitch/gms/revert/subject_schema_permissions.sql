-- Revert subject_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA subject FROM "gms-user";

COMMIT;
