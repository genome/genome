-- Verify model.build_input.index_build_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'm_bi_build_id_index';

ROLLBACK;
