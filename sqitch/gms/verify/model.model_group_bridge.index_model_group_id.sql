-- Verify model.model_group_bridge.index_model_group_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_group_bridge_model_group_id';

ROLLBACK;
