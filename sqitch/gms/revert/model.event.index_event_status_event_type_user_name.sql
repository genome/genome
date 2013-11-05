-- Revert model.event.index_event_status_event_type_user_name

BEGIN;

DROP INDEX model.idx_m_e_es_et_un;

COMMIT;
