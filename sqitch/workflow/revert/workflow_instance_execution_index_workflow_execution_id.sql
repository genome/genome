-- Revert workflow_instance_execution_index_workflow_execution_id

BEGIN;

DROP INDEX workflow.instance_execution_workflow_execution_id_idx;

COMMIT;
