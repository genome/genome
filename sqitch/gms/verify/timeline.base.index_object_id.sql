-- Verify timeline.base.index_object_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'base_object_id_idx';

ROLLBACK;
