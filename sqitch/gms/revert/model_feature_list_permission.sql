-- Revert model_feature_list_permission

BEGIN;

REVOKE ALL ON TABLE model.feature_list FROM "gms-user";

COMMIT;
