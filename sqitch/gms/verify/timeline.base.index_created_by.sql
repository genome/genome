-- Verify timeline.base.index_created_by

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'base_created_by_idx';

ROLLBACK;
