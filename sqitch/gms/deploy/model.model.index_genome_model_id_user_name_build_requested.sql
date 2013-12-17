-- Deploy model.model.genome_model_id_user_name_build_requested
-- requires: model_model

BEGIN;

CREATE INDEX m_m_id_user_name_build_requested_index on model.model using btree (genome_model_id, user_name, build_requested);

COMMIT;
