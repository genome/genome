-- Revert model_model_link_permission

BEGIN;

REVOKE ALL ON TABLE model.model_link FROM "gms-user";

COMMIT;
