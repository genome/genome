-- Revert result_param_permission

BEGIN;

REVOKE ALL ON TABLE result.param FROM "gms-user";

COMMIT;
