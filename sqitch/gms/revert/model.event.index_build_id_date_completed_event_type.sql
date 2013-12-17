-- Revert model.event.index_build_id_date_completed_event_type

BEGIN;

DROP INDEX model.idx_m_e_bi_dc_et;

COMMIT;
