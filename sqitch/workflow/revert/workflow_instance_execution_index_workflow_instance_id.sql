-- Revert workflow_instance_execution_index_workflow_instance_id

BEGIN;

DROP INDEX workflow.instance_execution_workflow_instance_id_idx;

COMMIT;
