-- Verify timeline.allocation.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'allocation_name_idx';

ROLLBACK;
