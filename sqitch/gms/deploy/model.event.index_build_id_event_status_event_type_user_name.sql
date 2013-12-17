-- Deploy model.event.build_id_event_status_event_type_user_name
-- requires: model_event

BEGIN;

CREATE INDEX idx_m_e_bi_es_et_un on model.event using btree (build_id, event_status, event_type, user_name);

COMMIT;
