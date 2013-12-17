-- Revert workflow_instance_index_name

BEGIN;

DROP INDEX workflow.instance_name_idx;

COMMIT;
