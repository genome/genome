-- Revert model.model_group_bridge.index_model_group_id

BEGIN;

DROP INDEX model.model_group_bridge_model_group_id;

COMMIT;
