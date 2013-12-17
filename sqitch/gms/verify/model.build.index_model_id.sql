-- Verify model.build.index_model_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'build_model_index';

ROLLBACK;
