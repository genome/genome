-- Verify workflow_instance_index_parent_instance_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instance_parent_instance_id_idx';

ROLLBACK;
