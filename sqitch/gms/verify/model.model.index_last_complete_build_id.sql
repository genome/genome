-- Verify model.model.index_last_complete_build_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'm_m_last_complete_build_id_index';

ROLLBACK;
