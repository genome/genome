-- Revert instrument_schema_permissions

BEGIN;

REVOKE ALL ON SCHEMA instrument FROM "gms-user";

COMMIT;
