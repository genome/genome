-- Deploy process_process_permissions
-- requires: process_process

BEGIN;

REVOKE DELETE ON TABLE process.process FROM "gms-user";

COMMIT;
