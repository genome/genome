-- Revert result.param.index_param_name

BEGIN;

DROP INDEX result.srp_pn2;

COMMIT;
