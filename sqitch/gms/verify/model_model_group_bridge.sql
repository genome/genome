-- Verify model_model_group_bridge

BEGIN;

SELECT model_group_id, model_id
FROM model.model_group_bridge
WHERE FALSE;

ROLLBACK;
