-- Deploy model.build_input.build_id
-- requires: model_build_input

BEGIN;

CREATE INDEX m_bi_build_id_index on model.build_input using btree (build_id);

COMMIT;
