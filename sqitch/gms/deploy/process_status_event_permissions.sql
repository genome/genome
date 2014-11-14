-- Deploy process_status_event_permissions
-- requires: process_status_event

BEGIN;

REVOKE DELETE,UPDATE ON TABLE process.status_event FROM "gms-user";

COMMIT;
