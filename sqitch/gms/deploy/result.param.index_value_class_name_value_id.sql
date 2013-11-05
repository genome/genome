-- Deploy result.param.value_class_name_value_id
-- requires: result_param

BEGIN;

CREATE INDEX result_param_value_class_id_index on result.param using btree (value_class_name, value_id);

COMMIT;
