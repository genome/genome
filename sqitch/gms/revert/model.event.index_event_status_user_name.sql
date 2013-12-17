-- Revert model.event.index_event_status_user_name

BEGIN;

DROP INDEX model.idx_m_e_es_un;

COMMIT;
