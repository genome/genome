-- Revert result.param.index_value_id

BEGIN;

DROP INDEX result.result_param_value_id_index;

COMMIT;
