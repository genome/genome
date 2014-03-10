-- Revert model_event_permission

BEGIN;

REVOKE ALL ON TABLE model.event FROM "gms-user";

COMMIT;
