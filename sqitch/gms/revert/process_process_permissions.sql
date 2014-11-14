-- Revert process_process_permissions

BEGIN;

GRANT DELETE ON TABLE process.process TO "gms-user";

COMMIT;
