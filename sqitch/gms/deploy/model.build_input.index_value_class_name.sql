-- Deploy model.build_input.value_class_name
-- requires: model_build_input

BEGIN;

CREATE INDEX m_bi_value_class_name_index on model.build_input using btree (value_class_name);

COMMIT;
