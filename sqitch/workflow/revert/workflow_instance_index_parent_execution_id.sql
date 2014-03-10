-- Revert workflow_instance_index_parent_execution_id

BEGIN;

DROP INDEX workflow.instance_parent_execution_id_idx;

COMMIT;
