-- Verify model.model.index_genome_model_id_user_name_build_requested

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'm_m_id_user_name_build_requested_index';

ROLLBACK;
