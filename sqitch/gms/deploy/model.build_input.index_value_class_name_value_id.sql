-- Deploy model.build_input.value_class_name_value_id
-- requires: model_build_input

BEGIN;

CREATE INDEX idx_m_bi_vcn_vi on model.build_input using btree (value_class_name, value_id);

COMMIT;
