-- Revert model_build_permission

BEGIN;

REVOKE ALL ON TABLE model.build FROM "gms-user";

COMMIT;
