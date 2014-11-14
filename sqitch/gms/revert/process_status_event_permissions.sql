-- Revert process_status_event_permissions

BEGIN;

GRANT DELETE,UPDATE ON TABLE process.status_event TO "gms-user";

COMMIT;
