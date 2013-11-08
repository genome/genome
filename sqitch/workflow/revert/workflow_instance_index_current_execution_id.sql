-- Revert workflow_instance_index_current_execution_id

BEGIN;

DROP INDEX workflow.instance_current_execution_id_idx;

COMMIT;
