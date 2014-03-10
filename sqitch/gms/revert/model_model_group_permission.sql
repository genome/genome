-- Revert model_model_group_permission

BEGIN;

REVOKE ALL ON TABLE model.model_group FROM "gms-user";

COMMIT;
