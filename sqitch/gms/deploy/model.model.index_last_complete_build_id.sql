-- Deploy model.model.last_complete_build_id
-- requires: model_model

BEGIN;

CREATE INDEX m_m_last_complete_build_id_index on model.model using btree (last_complete_build_id);

COMMIT;
