-- Revert model.model.index_current_running_build_id

BEGIN;

DROP INDEX model.m_m_current_running_build_id_index;

COMMIT;
