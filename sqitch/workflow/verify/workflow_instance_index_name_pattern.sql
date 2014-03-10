-- Verify workflow_instance_index_name_pattern

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'workflow_instance_name_varchar_pattern_ops_index';

ROLLBACK;
