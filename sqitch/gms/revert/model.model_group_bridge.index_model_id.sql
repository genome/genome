-- Revert model.model_group_bridge.index_model_id

BEGIN;

DROP INDEX model.model_group_bridge_model_id;

COMMIT;
