-- Revert model.event.index_parent_event_id_ref_seq_id

BEGIN;

DROP INDEX model.idx_m_e_pei_rsi;

COMMIT;
