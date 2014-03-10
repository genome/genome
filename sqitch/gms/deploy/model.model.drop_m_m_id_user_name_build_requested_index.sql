-- Deploy model.model.drop_m_m_id_user_name_build_requested_index
-- requires: model.model.index_genome_model_id_user_name_build_requested

BEGIN;

DROP INDEX model.m_m_id_user_name_build_requested_index;

COMMIT;
