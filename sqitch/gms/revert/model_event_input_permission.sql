-- Revert model_event_input_permission

BEGIN;

REVOKE ALL ON TABLE model.event_input FROM "gms-user";

COMMIT;
