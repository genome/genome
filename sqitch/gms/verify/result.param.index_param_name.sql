-- Verify result.param.index_param_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'srp_pn2';

ROLLBACK;
