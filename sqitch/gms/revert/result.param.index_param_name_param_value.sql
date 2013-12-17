-- Revert result.param.index_param_name_param_value

BEGIN;

DROP INDEX result.srp_pnpv2;

COMMIT;
