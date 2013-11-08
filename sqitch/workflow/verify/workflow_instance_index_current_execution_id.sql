-- Verify workflow_instance_index_current_execution_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instance_current_execution_id_idx';

ROLLBACK;
