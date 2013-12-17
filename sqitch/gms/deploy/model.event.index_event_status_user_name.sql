-- Deploy model.event.event_status_user_name
-- requires: model_event

BEGIN;

CREATE INDEX idx_m_e_es_un on model.event using btree (event_status, user_name);

COMMIT;
