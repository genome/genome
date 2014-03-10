-- Revert model.event.index_event_status_event_type_model_id

BEGIN;

DROP INDEX model.idx_m_e_es_et_mi;

COMMIT;
