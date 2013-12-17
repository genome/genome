-- Revert workflow_instance_index_name_pattern

BEGIN;

DROP INDEX workflow.workflow_instance_name_varchar_pattern_ops_index;

COMMIT;
