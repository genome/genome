-- Revert model_build_metric_permission

BEGIN;

REVOKE ALL ON TABLE model.build_metric FROM "gms-user";

COMMIT;
