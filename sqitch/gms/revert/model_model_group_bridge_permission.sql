-- Revert model_model_group_bridge_permission

BEGIN;

REVOKE ALL ON TABLE model.model_group_bridge FROM "gms-user";

COMMIT;
