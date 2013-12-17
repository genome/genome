-- Revert model.model.index_build_requested

BEGIN;

DROP INDEX model.m_m_build_requested_index;

COMMIT;
