-- Revert model.model.index_last_complete_build_id

BEGIN;

DROP INDEX model.m_m_last_complete_build_id_index;

COMMIT;
