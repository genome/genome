-- Verify workflow_instance_index_workflow_plan_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instance_workflow_plan_id_idx';

ROLLBACK;
