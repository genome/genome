-- Revert result.param.index_software_result_id

BEGIN;

DROP INDEX result.srp_sri2;

COMMIT;
