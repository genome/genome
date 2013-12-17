-- Revert result_metric_permission

BEGIN;

REVOKE ALL ON TABLE result.metric FROM "gms-user";

COMMIT;
