-- Verify workflow_execution_metric_index_workflow_execution_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'execution_metric_workflow_execution_id_idx';

ROLLBACK;
