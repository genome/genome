-- Revert model.model.drop_m_m_id_user_name_build_requested_index

BEGIN;

CREATE INDEX m_m_id_user_name_build_requested_index on model.model using btree (genome_model_id, user_name, build_requested);

COMMIT;
