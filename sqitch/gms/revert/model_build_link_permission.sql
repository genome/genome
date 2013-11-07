-- Revert model_build_link_permission

BEGIN;

REVOKE ALL ON TABLE model.build_link FROM "gms-user";

COMMIT;
