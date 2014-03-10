-- Revert result.param.index_name

BEGIN;

DROP INDEX result.result_param_name_index;

COMMIT;
