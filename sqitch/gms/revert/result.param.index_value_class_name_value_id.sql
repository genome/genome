-- Revert result.param.index_value_class_name_value_id

BEGIN;

DROP INDEX result.result_param_value_class_id_index;

COMMIT;
