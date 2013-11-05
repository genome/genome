-- Verify result.param.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'result_param_name_index';

ROLLBACK;
