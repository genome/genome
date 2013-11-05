-- Verify result.param.index_software_result_id_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'result_param_id_name';

ROLLBACK;
