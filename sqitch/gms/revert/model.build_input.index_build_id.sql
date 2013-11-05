-- Revert model.build_input.index_build_id

BEGIN;

DROP INDEX model.m_bi_build_id_index;

COMMIT;
