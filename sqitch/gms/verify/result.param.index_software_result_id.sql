-- Verify result.param.index_software_result_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'srp_sri2';

ROLLBACK;
