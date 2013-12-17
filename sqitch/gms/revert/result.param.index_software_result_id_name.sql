-- Revert result.param.index_software_result_id_name

BEGIN;

DROP INDEX result.result_param_id_name;

COMMIT;
