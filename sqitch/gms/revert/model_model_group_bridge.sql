-- Revert model_model_group_bridge

BEGIN;

DROP TABLE IF EXISTS model.model_group_bridge;

COMMIT;
