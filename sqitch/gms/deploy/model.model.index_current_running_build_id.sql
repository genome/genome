-- Deploy model.model.current_running_build_id
-- requires: model_model

BEGIN;

CREATE INDEX m_m_current_running_build_id_index on model.model using btree (current_running_build_id);

COMMIT;
