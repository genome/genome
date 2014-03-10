-- Revert model_model_input_permission

BEGIN;

REVOKE ALL ON TABLE model.model_input FROM "gms-user";

COMMIT;
