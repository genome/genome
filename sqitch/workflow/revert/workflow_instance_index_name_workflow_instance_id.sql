-- Revert workflow_instance_index_name_workflow_instance_id

BEGIN;

DROP INDEX workflow.instance_name_workflow_instance_id_idx;

COMMIT;
