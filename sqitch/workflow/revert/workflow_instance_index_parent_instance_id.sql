-- Revert workflow_instance_index_parent_instance_id

BEGIN;

DROP INDEX workflow.instance_parent_instance_id_idx;

COMMIT;
