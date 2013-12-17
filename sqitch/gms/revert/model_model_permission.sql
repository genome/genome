-- Revert model_model_permission

BEGIN;

REVOKE ALL ON TABLE model.model FROM "gms-user";

COMMIT;
