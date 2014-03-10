-- Verify workflow_instance_index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instance_name_idx';

ROLLBACK;
