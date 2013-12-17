-- Revert workflow_schema_permission

BEGIN;

REVOKE ALL ON SCHEMA workflow FROM "gms-user";

COMMIT;
