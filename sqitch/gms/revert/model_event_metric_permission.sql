-- Revert model_event_metric_permission

BEGIN;

REVOKE ALL ON TABLE model.event_metric FROM "gms-user";

COMMIT;
