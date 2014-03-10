-- Verify timeline.allocation.index_absolute_path

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'allocation_absolute_path_idx';

ROLLBACK;
