-- Verify result.param.index_param_name_param_value

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'srp_pnpv2';

ROLLBACK;
