-- Verify workflow_instance_execution_index_workflow_execution_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instance_execution_workflow_execution_id_idx';

ROLLBACK;
