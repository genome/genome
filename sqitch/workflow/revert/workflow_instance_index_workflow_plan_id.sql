-- Revert workflow_instance_index_workflow_plan_id

BEGIN;

DROP INDEX workflow.instance_workflow_plan_id_idx;

COMMIT;
