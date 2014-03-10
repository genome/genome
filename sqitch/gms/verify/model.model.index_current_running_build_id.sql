-- Verify model.model.index_current_running_build_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'm_m_current_running_build_id_index';

ROLLBACK;
