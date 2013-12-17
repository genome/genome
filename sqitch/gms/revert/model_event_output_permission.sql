-- Revert model_event_output_permission

BEGIN;

REVOKE ALL ON TABLE model.event_output FROM "gms-user";

COMMIT;
