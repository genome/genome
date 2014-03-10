-- Revert model_build_input_permission

BEGIN;

REVOKE ALL ON TABLE model.build_input FROM "gms-user";

COMMIT;
