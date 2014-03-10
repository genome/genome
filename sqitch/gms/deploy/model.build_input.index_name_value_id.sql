-- Deploy model.build_input.name_value_id
-- requires: model_build_input

BEGIN;

CREATE INDEX idx_m_bi_n_vi on model.build_input using btree (name, value_id);

COMMIT;
